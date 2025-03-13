#!/usr/bin/env python
from larixite.webapp import app

app.jinja_env.cache = {}
app.run(debug=True, port=11564)
