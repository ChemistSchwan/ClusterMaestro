coverage run --source ClusterMaster -m pytest
coverage report --skip-covered --omit=*setup.py -m
