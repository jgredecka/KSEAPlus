import io
from io import StringIO, BytesIO
import os
import sys
from datetime import timedelta
import importlib
import pandas as pd
import urllib.parse

import redis
from flask import Flask, flash, render_template, url_for, session, redirect, request, send_file, send_from_directory, jsonify
from flask_kvsession import KVSessionExtension
from simplekv.memory.redisstore import RedisStore

from flask_wtf import FlaskForm
from wtforms import SelectField, IntegerField

from databases import uploadPSP, uploadPDTS, uploadEDGES
import appconfig
from celery import Celery
from celery.result import AsyncResult
from celery.signals import task_postrun

# Set of allowed file extensions.
ALLOWED_EXTENSIONS = set(['tsv'])

sys.path.append(os.getcwd())

# Create a redis instance with connection pool.
redis_url = urllib.parse.urlparse(appconfig.REDIS_URL)
store = RedisStore(redis.StrictRedis(host = redis_url.hostname, port = redis_url.port, password = redis_url.password, db = 0, socket_connect_timeout = 30, socket_timeout = 30, socket_keepalive=True, retry_on_timeout=True, health_check_interval=55))
app = Flask(__name__)
KVSessionExtension(store, app)

# App configuration settings.
app.permanent_session_lifetime = appconfig.PERMANENT_SESSION_LIFETIME
app.secret_key = appconfig.SECRET_KEY
app.debug = appconfig.DEBUG

# Celery object for background tasks.
def make_celery(app):
    celery = Celery(
        app.import_name,
        backend=appconfig.CELERY_RESULT_BACKEND,
        broker=appconfig.CELERY_BROKER_URL
    )
    celery.conf.update(app.config)
    
    class ContextTask(celery.Task):
        def __call__(self, *args, **kwargs):
            with app.app_context():
                return self.run(*args, **kwargs)

    celery.Task = ContextTask
    return celery

app.config.update(
    CELERY_BROKER_URL = appconfig.REDIS_URL,
    CELERY_RESULT_BACKEND = appconfig.REDIS_URL,
    CELERY_IGNORE_RESULT = False,
    CELERY_TASK_TRACK_STARTED = True,
    CELERY_TASK_ALWAYS_EAGER = False,
    CELERY_TASK_EAGER_PROPAGATES = False,
    CELERY_TASK_SERIALIZER = 'json',
    CELERY_ACCEPT_CONTENT = ['json'],
    BROKER_POOL_LIMIT = 1,
    CELERYD_MAX_TASKS_PER_CHILD = 1,
    BROKER_CONNECTION_MAX_RETRIES = None
)
celery = make_celery(app)

# Create a celery task for ksea analyses to run in background.
@celery.task(name='tasks.run.runAlg')
def runAlg(script, ks_db, graphics, df, min_sub):
    mod = importlib.import_module(script)
    return mod.userInput(ks_db, graphics, df, min_sub)

# Attempts to free memory after each task.
@task_postrun.connect 
def gc_after_task(**kwargs):
    import gc
    gc.collect() 

# Upload the PhosphoSitePlus database.
psp_db = uploadPSP()
# Upload the PDTS database.
pdts_db = uploadPDTS()
# Upload the EDGES database.
edges_db = uploadEDGES()

# A dictionary mapping the database name to the array storing the database entries.
# Used for the conditional loop on the 'upload' page.
db_map = {"psp": psp_db, "pdts": pdts_db, "edges": edges_db}
single_list = ["ztest_single", "karp_single", "ks_single"]
multi_list = ["ztest_multi", "karp_multi", "ks_multi"]

# Function that checks if a file extension is valid. Must return the boolean 'true' to proceed.
def allowed_file(filename):
    return '.' in filename and filename.rsplit('.', 1)[1].lower() in ALLOWED_EXTENSIONS

# Selection form for available databases is defined.
class DbForm(FlaskForm):
    choices = [("psp", "PhosphoSitePlus"), ("pdts", "PDTS"), ("edges", "EDGES")]
    select_db = SelectField(choices=choices)
    
# Selection form for all KSEA algorithms is defined here.
class AlgForm(FlaskForm):
    choices = [("ks_single", "Kolmogorov-Smirnov (Single Sample)"), ("ks_multi", "Kolmogorov-Smirnov (Multiple Samples)"), ("karp_single", "KARP (Single Sample)"), ("karp_multi", "KARP (Multiple Samples)"), ("ztest_single", "Z-test (Single Sample)"), ("ztest_multi", "Z-test (Multiple Samples)")]
    select_alg = SelectField(choices=choices)
    
# Selection form for whether or not to plot a graph.
class GraphicsForm(FlaskForm):
    choices = [("no", "No"), ("yes", "Yes")]
    select_graphics = SelectField(choices=choices)

# Plot parameters: minimum number of substrates per kinase.
class SubForm(FlaskForm):
    sub_choice = IntegerField("substrate_count", default=5)

# Home page route.
@app.route("/")
def index():
    return render_template("index.html", title="Home")
    
# This route allows the user to upload a .tsv file and set various parameters.
# A celery task is created and runAlg runs the algorithm in the bg.
@app.route("/upload", methods = ['GET', 'POST'])
def upload():
    db_form = DbForm()
    alg_form = AlgForm()
    sub_form = SubForm()
    plot_form = GraphicsForm()
    select_db = None
    select_alg = None
    sub_choice = None
    select_graphics = None
    if request.method == 'POST':
        if int(request.content_length) > 5 * 1024 * 1024:
            flash('Files larger than 5MB are not allowed.')
            return redirect(url_for('upload'))
        file = request.files['file']
        if file and allowed_file(file.filename):
            session.permanent = True
            session.regenerate()
            f = file.read()
            stream = BytesIO()
            stream.write(f)
            stream.seek(0)
            df = pd.read_csv(stream, sep="\t")
            df=df.to_json(orient="split")
            stream.close()
            select_db = db_form.select_db.data
            select_alg = alg_form.select_alg.data
            min_sub = sub_form.sub_choice.data
            graphics = plot_form.select_graphics.data
            alg_type = select_alg.split("_")[1]
            if alg_type == "single":
                alg_list = single_list
                for x in alg_list:
                    if select_alg == x:
                        script = x
                        ks_db = db_map[select_db]
                        res = runAlg.delay(script, ks_db, graphics, df, min_sub)
                        taskid = res.task_id
                        return redirect(url_for('show_results', alg_type=alg_type, taskid=taskid))
                    else:
                        continue
            elif alg_type == "multi":
                alg_list = multi_list
                for x in alg_list:
                    if select_alg == x:
                        script = x
                        ks_db = db_map[select_db]
                        res = runAlg.delay(script, ks_db, graphics, df, min_sub)
                        taskid = res.task_id
                        return redirect(url_for('show_results', alg_type=alg_type, taskid=taskid))
                    else:
                        continue
    return render_template("upload.html", title="Upload File", db_form=db_form, alg_form=alg_form, sub_form=sub_form, plot_form=plot_form)

@app.route("/getting-started")   
def getting_started():
    return render_template("getting_started.html", title="Getting Started")
    
@app.route("/algorithms")
def algorithms():
    return render_template("algorithms.html", title="Algorithm Overview")

@app.route("/contact")
def contact():
    return render_template("contact.html", title="Contact")

@app.route("/sitemap.xml")
def sitemap():
    return send_from_directory('static', filename='sitemap/sitemap.xml')

@app.route("/test-results", methods=['GET', 'POST'])
def get_results():
    tries = 3
    for attempt in range(tries):
        try:
            alg = request.args.get('alg')
            taskid = request.args.get('taskid')
            r = celery.AsyncResult(taskid)
            if r.ready() == True and alg == 'single':
                try:
                    session.regenerate()
                    result = r.result
                    r.forget()
                    scores = result[0]
                    links = result[1]
                    plot = result[2]
                    session["scores"] = scores
                    session["links"] = links
                    session["plot"] = plot
                    sdf = pd.read_json(scores, orient='split')
                    ldf = pd.read_json(links, orient='split')
                    score_tab=sdf.to_html(index=False, classes = 'row-border hover stripe" id = "resTable').replace('border="1"','border="0"')
                    links_tab = ldf.to_html(index=False, classes = 'row-border hover stripe" id = "linkTable').replace('border="1"','border="0"')
                    return jsonify({'plot':plot, 'scores':score_tab, 'links':links_tab, 'status': 'complete'})
                except(ValueError, TypeError):
                    r.forget()
                    return jsonify({'err':'An error has occured. Please check your dataset and plot parameters.', 'status': 'error'})
            elif r.ready() == True and alg == 'multi':
                try:
                    session.regenerate()
                    result = r.result
                    r.forget()
                    scores = result[0]
                    links = result[1]
                    plot = result[2]
                    session["scores"] = scores
                    session["links"] = links
                    session["plot"] = plot
                    return jsonify({'plot':plot, 'scores': 'Kinase activity data available only as a downloadable .csv file.', 'links': 'K-S relationship data available only as a downloadable .csv file.', 'status': 'complete'})
                except(ValueError, TypeError):
                    r.forget()
                    return jsonify({'err':'An error has occured. Please check your dataset and plot parameters.', 'status': 'error'})
            else:
                return jsonify({'status':'pending'})
        except Exception as e:
            if attempt < tries - 1:
                continue
            else:
                raise

@app.route("/ksea/<alg_type>/<taskid>", methods=['GET', 'POST'])
def show_results(alg_type, taskid):
    uid = taskid
    placeholder=''
    return render_template('test.html', placeholder=placeholder, alg_type=alg_type, taskid=taskid, title="KSEA Results", uid=uid)

# send_file requires BytesIO.
# Here SVG data is encoded to bytes and written into the buffer.
@app.route("/download/fig/<uid>")
def fig_download(uid):
    try:
        graph = session.get("plot")
        buffer = BytesIO()
        buffer.write(graph.encode('utf-8'))
        buffer.seek(0)
        return send_file(buffer, mimetype='image/svg+xml', attachment_filename="plot-"+uid+".svg", as_attachment=True)
    except AttributeError:
        return render_template("timeout.html", title="Session expired")

# Here kinase-score json string is converted back into a dataframe and then written into a proxy string buffer as a csv file. It is then encoded to downloadable bytes data.
@app.route("/download/scores/<uid>")
def download_scores(uid):
    try:
        scores = session.get("scores")
        scores = pd.read_json(scores, orient='split')
        proxy = StringIO()
        scores.to_csv(proxy, index=False)
        proxy.seek(0)
        data=proxy.getvalue()
        proxy.close()
        byte_data=data.encode('utf-8')
        buffer = BytesIO()
        buffer.write(byte_data)
        buffer.seek(0)
        return send_file(buffer, mimetype='text/csv', attachment_filename="ksea_scores-"+uid+".csv", as_attachment=True)
    except ValueError:
        return render_template("timeout.html", title="Session expired")
    
# Here, Kinase-Substrate-Relationships data is manipulated as above in order to be downloadable.
@app.route("/download/links/<uid>")
def download_links(uid):
    try:
        links = session.get("links")
        links = pd.read_json(links, orient='split')
        proxy = StringIO()
        links.to_csv(proxy, index=False)
        proxy.seek(0)
        data = proxy.getvalue()
        proxy.close()
        byte_data = data.encode('utf-8')
        buffer = BytesIO()
        buffer.write(byte_data)
        buffer.seek(0)
        return send_file(buffer, mimetype='text/csv', attachment_filename="ks-links-"+uid+".csv", as_attachment=True)
    except ValueError:
        return render_template("timeout.html", title="Session expired")
    
# Run in production
if __name__ == '__main__':
    app.run()