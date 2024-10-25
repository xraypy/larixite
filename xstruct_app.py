#!/usr/bin/env python
from cif4xas import app

app.jinja_env.cache = {}
app.run(debug=True, port=11564)
