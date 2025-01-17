# KSEAPlus

KSEAPlus is a computational tool that implements a number of kinase activity inference algorithms to allow for the quantification of kinase activity from mass spectrometry-based phosphoproteomic data. KSEAPlus enables the analysis of multi-condition datasets. Currently offered algorithms include KSEA Z-Test, KARP and Kolmogorov-Smirnov test. All algorithms can be ran against one of three kinase-substrate databases: PhosphoSitePlus, PDTS and EDGES.

## Availability

KSEAPlus is currently running on a free Heroku dyno: https://ksea-plus.herokuapp.com

## Built With
* Flask: serving as a Python web microframework.
* Celery: background task manager for running the analyses in application background.
* Gunicorn: acting as a WSGI server.

## Documentation
Full documentation containing instructions on how to use KSEAPlus is available [here](static/starter-pack/KSEAPlus-User-Guide.pdf).
