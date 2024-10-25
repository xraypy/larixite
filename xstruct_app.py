#!/usr/bin/env python
from xstructures.webapp import app

app.jinja_env.cache = {}
app.run(debug=True, port=11564)
