#!/usr/bin/env python
import os
import time
import numpy as np
from pathlib import Path

from flask import (Flask, redirect, url_for, render_template, flash,
                   request, session, Response, send_from_directory)
from werkzeug.utils import secure_filename

from larixite import get_amcsd, cif_cluster, cif2feffinp, read_cif_structure
from larixite.utils import get_homedir

from xraydb import atomic_number
from xraydb.chemparser import chemparse

top =  Path(__file__).absolute().parent

app = Flask('larixite',
            static_folder=Path(top, 'static'),
            template_folder=Path(top, 'templates'))


UPLOAD_FOLDER = get_homedir()
ALLOWED_EXTENSIONS = {'cif', }

app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER
app.config['MAX_CONTENT_LENGTH'] = 2**24
app.secret_key = 'persuasive butternut lanterns'

app.config.from_object(__name__)

cifdb = config = None

def connect(session, clear=False, cifid=None):
    global cifdb, config
    if cifdb is None:
        cifdb = get_amcsd()

    if config is None or clear:
        config = {'cifdict': {},
                  'with_h': 0,
                  'edges': ['K', 'L3', 'L2', 'L1', 'M5', 'M4', 'M3'],
                  'cluster_size': 7.0,
                  'mineral': '',
                  'elems_in': '',
                  'elems_out': '',
                  'strict': True,
                  'absorber': '',
                  'edge': '',
                  'ciftext': ''}

    if cifid is None:
        mode = 'browse'
        config['ciftext'] = ''
        config['ciffile'] = ''
    else:
        if cifid.endswith('.cif') or cifid.endswith('.txt'):
            mode = 'file'
        else:
            try:
                cifid = int(cifid)
                mode = 'amsid'
            except ValueError:
                mode = 'file'

    if mode == 'amsid':
        config['cifid'] = int(cifid)
        cif = cifdb.get_cif(cifid)
        config['ciftext'] = cif.ciftext.strip()
        config['ciffile'] = ''
    elif mode == 'file':
        config['cifid'] = cifid
        config['ciffile'] = cifid
        try:
            path =Path(app.config["UPLOAD_FOLDER"], cifid).absolute()
            if path.exists():
                with open(path, 'r') as fh:
                    config['ciftext'] = fh.read().strip()
        except:
            pass

def allowed_file(filename):
    return '.' in filename and \
           filename.rsplit('.', 1)[1].lower() in ALLOWED_EXTENSIONS


@app.route('/favicon.ico')
def favicon():
    return send_from_directory(app.static_folder, 'ixas_logo.ico',
                               mimetype='image/vnd.microsoft.icon')

@app.route('/cifs')
@app.route('/cifs', methods=['GET', 'POST'])
@app.route('/cifs/<cifid>', methods=['GET', 'POST'])
def cifs(cifid=None):
    global cifdb, config
    connect(session, clear=False, cifid=cifid)

    if len(config['ciftext']) > 4:
        config['ciftext'] = '\n\n' + config['ciftext']

        try:
            t = read_cif_structure(config['ciftext'])
        except ValueError:
            error = f"could not read CIF:   {conf['cifid']}"
            return render_template('index.html',
                            error=error)
        cluster = cif_cluster(config['ciftext'])
        config['atoms'] = []
        config['cif_link'] = None
        config['feff_links'] = {}
        config['all_sites'] = cluster.all_sites

        for at in chemparse(cluster.formula.replace(' ', '')):
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
                    pass
            config['cifdict'] = cifdict
            config['cifid'] = cifid  = None
            config['ciftext'] = ''
            config['ciffile']  = ''
            config['atoms'] = []
            config['cif_link'] = None
            config['feff_links'] = {}


        elif 'feff' in request.form.keys():
            config['absorber'] = absorber = request.form.get('absorbing_atom')
            config['edge'] = edge =request.form.get('edge')
            with_h = request.form.get('with_h') in ('1', 'on', 'True', 1, True)
            config['with_h'] = with_h = int(with_h)
            config['cluster_size'] = cluster_size = request.form.get('cluster_size')

            if config['ciftext'] is not None:
                cluster = cif_cluster(config['ciftext'], absorber=absorber)
                config['all_sites'] = cluster.all_sites
                cifid = config['cifid']
                if config['ciffile'] not in ('', None, 'None'):
                    cifid = config['ciffile']
                elif config['cifid'] is not None:
                    cifid = config['cifid']
                    config['cif_link'] = (f'AMSCD_{cifid}.cif', cifid)
                config['feff_links'] = {}
                for slabel, sindex in config['all_sites'][absorber].items():
                    link = f'feff_{cifid}_{absorber}{sindex}_{edge}.inp'
                    config['feff_links'][link] = (slabel, cifid, absorber, edge, sindex,
                                                  cluster_size, with_h)

    return render_template('index.html', **config)


@app.route('/')
def index():
    global cifdb, config
    connect(session, clear=True)
    return render_template('index.html', **config)


@app.route('/feffinp/<cifid>/<absorber>/<site>/<edge>/<cluster_size>/<with_h>/<fname>')
def feffinp(cifid=None, absorber=None, site=1, edge='K', cluster_size=7.0,
            with_h=False, fname=None):
    connect(session, cifid=cifid)
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
        if ( ('M' in absorber)
             and (not absorber.startswith('Mo'))
             and (not absorber.startswith('Mn'))):
            try:
                stoich = chemparse(absorber[:].replace('M', ''))
            except:
                stoich = {}

    if len(stoich) != 1:
        txt = f"could not parse absorber: '{absorber}' for CIF {cifid}  {stoich}"
    else:
        try:
            xcifid = int(cifid)
        except:
            xcifid = None

        absorber = list(stoich.keys())[0]
        if with_h in  ('1', 'on', 'True', 1, True):
            with_h = 1
        else:
            with_h = 0
        feffinp = cif2feffinp(config['ciftext'], absorber, edge=edge,
                    cluster_size=float(cluster_size), cifid=xcifid,
                    with_h=with_h, absorber_site=int(site))

    return Response(feffinp, mimetype='text/plain')


@app.route('/ciffile/<cifid>/<fname>')
def ciffile(cifid=None, fname='amcsd.cif'):
    global cifdb, config
    connect(session, cifid=cifid)
    return Response(config['ciftext'], mimetype='text/plain')

@app.route('/upload/')
def upload():
    return render_template('upload.html')

@app.route('/upload_cif', methods=['GET', 'POST'])
def upload_cif():
    if request.method == 'POST':
        if 'file' not in request.files:
            flash('No file part')
            return redirect(request.url)
        file = request.files['file']
        if file.filename == '':
            flash('No selected file')
            return redirect(request.url)
        if file and allowed_file(file.filename):
            filename = secure_filename(file.filename)
            path = Path(app.config['UPLOAD_FOLDER'], filename).absolute()
            file.save(path)
            time.sleep(0.25)
            error = None
            if path.exists():
                with open(path, 'r') as fh:
                    try:
                        t = read_cif_structure(fh.read())
                    except ValueError:
                        error = f"Pymatgen could not read CIF file: {filename}"
            else:
                error = f"could not upload CIF file: {filename}"
            if error is not None:
                return render_template('upload.html', error=error)

            return redirect(url_for('cifs', cifid=filename))

    return render_template('upload.html')



@app.route('/about/')
def about():
    return render_template('about.html')
