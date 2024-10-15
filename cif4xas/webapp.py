#!/usr/bin/env python
import os
import time
import numpy as np
from pathlib import Path

from flask import (Flask, redirect, url_for, render_template,
                   request, session, Response, send_from_directory)

import larch
from larch.xrd import get_amcsd
from larch.xrd.cif2feff import cif2feffinp, cif_sites
from xraydb import atomic_number
from xraydb.chemparser import chemparse

top =  Path(__file__).absolute().parent

app = Flask('cif4xas',
            static_folder=Path(top, 'static'),
            template_folder=Path(top, 'templates'))

app.config.from_object(__name__)

cifdb = config = None

def connect(session, force_refresh=False):
    global cifdb, config
    if cifdb is None:
        cifdb = get_amcsd()

    if config is None:
        config = {'cifid': None,
                   'mineral': '',
                   'elems_in': '',
                   'elems_out': '',
                   'strict': True,
                   'cifdict': {},
                   'atoms': [],
                   'with_h': False,
                   'edges': ['K', 'L3', 'L2', 'L1', 'M5', 'M4', 'M3'],
                   'cluster_size': 7.0,
                   'ciftext' : '',
                   'absorber': '',
                   'edge': '',
                   'feff_fname': '',
                   'fdmnes_fname': '',
        }


@app.route('/favicon.ico')
def favicon():
    return send_from_directory(app.static_folder, 'ixas_logo.ico',
                               mimetype='image/vnd.microsoft.icon')

@app.route('/')
@app.route('/', methods=['GET', 'POST'])
@app.route('/<cifid>/', methods=['GET', 'POST'])
def index(cifid=None):
    connect(session)
    global cifdb, config
    if cifid is not None:
        config['cifid'] = cifid

        cif = cifdb.get_cif(cifid)
        config['ciftext'] = f"\n{cif.ciftext.strip()}"
        config['atoms'] = []
        config['feff_fname'] = ''
        for at in cif.atoms_sites:
            if at not in ('H', 'D') and at not in config['atoms']:
                config['atoms'].append(at)

    if request.method == 'POST':
        if 'search' in request.form.keys():
            mineral = request.form.get('mineral')
            elems_in = request.form.get('elems_in').strip()
            elems_out = request.form.get('elems_out').strip()
            strict = request.form.get('strict') is not None
            config.update({'elems_in': elems_in, 'elems_out': elems_out,
                           'mineral': mineral, 'strict': strict})
            if len(mineral) > 2 and not mineral.endswith('*'):
                mineral = mineral + '*'

            contains_elements = excludes_elements = None
            if len(elems_in) >  0:
                contains_elements = [a.strip().title() for a in elems_in.split(',')]
            if len(elems_out) >  0:
                excludes_elements = [a.strip().title() for a in elems_out.split(',')]

            cifdict = {}
            all_cifs = cifdb.find_cifs(mineral_name=mineral,
                                       contains_elements=contains_elements,
                                       excludes_elements=excludes_elements,
                                       strict_contains=strict, max_matches=500)

            for cif in all_cifs:
                try:
                    label = cif.formula.replace(' ', '')
                    mineral = cif.get_mineralname()
                    year = cif.publication.year
                    journal= cif.publication.journalname
                    cid = cif.ams_id
                    label = f'{label}: {mineral}'
                    cite = f'{journal} {year}'
                    if len(label) > 80 and len(cifict) > 125:
                        continue
                    cifdict[cid] = (label, cite)
                except:
                    print("could not add cif ", cif)
            config['cifdict'] = cifdict
            config['cifid'] = cifid  = None
            config['ciftext'] = ''
            config['atoms'] = []
            config['feff_fname'] = ''

        elif 'feff' in request.form.keys():
            config['absorber'] =  request.form.get('absorbing_atom')
            config['edge'] = request.form.get('edge')
            config['with_h'] = 1 if request.form.get('with_h') else 0
            config['cluster_size'] = request.form.get('cluster_size')

            absorber = config['absorber']
            edge = config['edge']
            cifid = config['cifid']
            config['feff_fname'] = f'feff_CIF{cifid}_{absorber}_{edge}.inp'

    return render_template('index.html', **config)


@app.route('/feffinp/<cifid>/<absorber>/<edge>/<cluster_size>/<with_h>/<fname>')
def feffinp(cifid=None, absorber=None, edge='K', cluster_size=7.0,
            with_h=False, fname=None):
    connect(session)
    global cifdb, config
    if absorber.startswith('Wat'):
        absorber.replace('Wat', 'O')
    if absorber.startswith('O-H'):
        absorber.replace('O-H', 'O')
    if absorber.startswith('D') and not absorber.startswith('Dy'):
        absorber.replace('D', 'H')


    try:
        stoich = chemparse(absorber)
    except ValueError:
s        if ( ('M' in absorber)
             and (not absorber.startswith('Mo'))
             and (not absorber.startswith('Mn'))):
            try:
                stoich = chemparse(absorber[:].replace('M', ''))
            except:
                stoich = {}

    if len(stoich) != 1:
        txt = f"could not parse absorber: '{absorber}' for CIF {cifid}  {stoich}"
    else:
        absorber = list(stoich.keys())[0]
        absorber_site = int(stoich[absorber]) + 1

        cif = cifdb.get_cif(cifid)

        title = f'*** feff input generated with cif4xas {time.ctime()}'
        feffinp = cif2feffinp(cif.ciftext, absorber, edge=edge,
                    cluster_size=float(cluster_size),
                    with_h=with_h, absorber_site=absorber_site)
        txt = title + '\n' + feffinp

    return Response(txt, mimetype='text/plain')

@app.route('/showcifuploaded/<fid>')
def showcifuploaded(fid=None):
    print('upload ')


@app.route('/upload/')
def upload():
    return render_template('upload.html')


@app.route('/about/')
def about():
    return render_template('about.html')
