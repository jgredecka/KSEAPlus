import io
from io import StringIO, BytesIO
import os
import uuid
from datetime import timedelta
import importlib
import pandas as pd

import redis
from flask import Flask, flash, render_template, url_for, session, redirect, request, send_file
from flask_kvsession import KVSessionExtension
from simplekv.memory.redisstore import RedisStore

from flask_wtf import FlaskForm
from wtforms import SelectField, IntegerField

from databases import uploadPSP, uploadPDTS, uploadEDGES
import appconfig

# Set of allowed file extensions.
ALLOWED_EXTENSIONS = set(['tsv'])

store = RedisStore(redis.StrictRedis(host = '0.0.0.0', port = 6379, db = 0))
app = Flask(__name__)
KVSessionExtension(store, app)

# App configuration settings.
app.max_content_length = appconfig.MAX_CONTENT_LENGTH
app.permanent_session_lifetime = appconfig.PERMANENT_SESSION_LIFETIME
app.secret_key = appconfig.SECRET_KEY
app.debug = appconfig.DEBUG

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
db_list = ["psp", "pdts", "edges"]

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
    
class SubForm(FlaskForm):
    sub_choice = IntegerField("substrate_count", default=5)

# Home page route.
@app.route("/")
def index():
    return render_template("index.html", title="Home")
    
# This route allows the user to upload a .tsv file and set various parameters.
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
        file = request.files['file']
        if file and allowed_file(file.filename):
            f = file.read()
            stream = BytesIO()
            stream.write(f)
            stream.seek(0)
            df = pd.read_csv(stream, sep="\t")
            stream.close()
            select_db = db_form.select_db.data
            select_alg = alg_form.select_alg.data
            min_sub = sub_form.sub_choice.data
            graphics = plot_form.select_graphics.data
            if select_alg.split("_")[1] == "single":
                alg_list = single_list
                for x in alg_list:
                    for y in db_list:
                        try:
                            if select_alg == x and select_db == y:
                                script = x + "_" + y
                                mod = importlib.import_module(script)
                                ks_db = db_map[y]
                                score_df, links_df, session["svg_plot"] = mod.userInput(df, min_sub, ks_db, graphics)
                                session["score_df"] = score_df.to_json(orient='split')
                                session["links_df"] = links_df.to_json(orient='split')
                                return redirect(url_for('single_results'))
                            else:
                                continue
                        except(ValueError, TypeError):
                            flash("An error has occured. Please check your dataset and plot parameters.")
                            return redirect(url_for('upload'))
            elif select_alg.split("_")[1] == "multi":
                alg_list = multi_list
                for x in alg_list:
                    for y in db_list:
                        try:
                            if select_alg == x and select_db == y:
                                script = x + "_" + y
                                mod = importlib.import_module(script)
                                ks_db = db_map[y]
                                score_df, links_df, session["svg_plot"] = mod.userInput(df, min_sub, ks_db, graphics)
                                session["score_df"] = score_df.to_json(orient='split')
                                session["links_df"] = links_df.to_json(orient='split')
                                return redirect(url_for('multi_results'))
                            else:
                                continue
                        except(ValueError, TypeError):
                            flash("An error has occured. Please check your dataset and plot parameters.")
                            return redirect(url_for('upload'))
    return render_template("upload.html", title="Upload File", db_form=db_form, alg_form=alg_form, sub_form=sub_form, plot_form=plot_form)

@app.route("/results/single", methods=['GET', 'POST'])
def single_results():
    barplot = session.get("svg_plot")
    score_data = session.get("score_df")
    score_data = pd.read_json(score_data, orient='split')
    links_data = session.get("links_df")
    links_data = pd.read_json(links_data, orient='split')
    return render_template("single_results.html", title="KSEA Results", score_data=score_data.to_html(index=False, classes = 'row-border hover stripe" id = "resTable').replace('border="1"','border="0"'), links_data=links_data.to_html(index=False, classes = 'row-border hover stripe" id = "linkTable').replace('border="1"','border="0"'), barplot=barplot)

@app.route("/results/multi")
def multi_results():
    heatmap = session.get("svg_plot")
    return render_template("multi_results.html", title="KSEA Results", heatmap=heatmap)

@app.route("/getting-started")   
def getting_started():
    return render_template("getting_started.html", title="Getting Started")
    
@app.route("/algorithms")
def algorithms():
    return render_template("algorithms.html", title="Algorithm Overview")

@app.route("/contact")
def contact():
    return render_template("contact.html", title="Contact")

# send_file requires BytesIO.
# Here SVG data is encoded to bytes and written into the buffer.
@app.route("/download/fig")
def fig_download():
    graph = session.get("svg_plot")
    buffer = BytesIO()
    buffer.write(graph.encode('utf-8'))
    buffer.seek(0)
    return send_file(buffer, mimetype='image/svg+xml', attachment_filename="plot.svg", as_attachment=True)

@app.route("/download/scores")
def download_scores():
    scores = session.get("score_df")
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
    return send_file(buffer, mimetype='text/csv', attachment_filename="ksea_scores.csv", as_attachment=True)
    
@app.route("/download/links")
def download_links():
    links = session.get("links_df")
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
    return send_file(buffer, mimetype='text/csv', attachment_filename="ks_links.csv", as_attachment=True)

if __name__ == '__main__':
    app.run(host='0.0.0.0')