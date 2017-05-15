from flask import Flask, render_template, request, jsonify, json, redirect, url_for, send_from_directory, Response
from werkzeug import secure_filename
import backend, readapi, readfiles, os
from wtforms import TextField, Form

# Initialize the Flask application

app = Flask(__name__)

class SearchForm(Form):
    autocomp = TextField('Insert Entrez_id', id='gene_autocomplete')

#backend.writeGeneList()
genes = backend.readGeneList()
genelist = ""
voilist = backend.createVoiList()
if not os.path.exists('uploads/'):
    os.makedirs('uploads/')
app.config['UPLOAD_FOLDER'] = 'uploads/'
app.config['ALLOWED_EXTENSIONS'] = set(['csv', 'gz'])

def allowed_file(filename):
    return '.' in filename and \
           filename.rsplit('.', 1)[1] in app.config['ALLOWED_EXTENSIONS']

@app.route('/')
def index():
    form = SearchForm(request.form)
    return render_template('index.html', form=form)

@app.route('/_jugex')
def jugex():
    #res = backend.performJugexFromFiles()
    res = backend.performJugex()
    return jsonify(result=res)

@app.route('/_getSelectedGenes')
def getSelectedGenes():
    global genelist
    gene = request.args.get('s')
    genelist += gene+','
    return jsonify(result=genelist)


@app.route('/uploadMultiple', methods=['POST'])
def uploadMultiple():
    form = SearchForm(request.form)
    if request.method == 'POST':
        files = request.files.getlist("file[]")
        for file in files:
            if file and allowed_file(file.filename):
                filename = secure_filename(file.filename)
                path = os.path.join(app.config['UPLOAD_FOLDER'], filename)
                file.save(path)
        return ('', 204)
        #return render_template('index.html', form=form)

@app.route('/_selectVois1', methods=['GET', 'POST'])
def selectVois1():
    form = SearchForm(request.form)
    if request.method == 'POST':
        region = request.form.get('regionList1')
        print(voilist[region])
        backend.voinames.append(voilist[region])
        return ('', 204)
        #return render_template('index.html', form=form)

@app.route('/_selectVois2', methods=['GET', 'POST'])
def selectVois2():
    form = SearchForm(request.form)
    if request.method == 'POST':
        region = request.form.get('regionList2')
        print(voilist[region])
        backend.voinames.append(voilist[region])
        return ('', 204)
        #return render_template('index.html', form=form)


@app.route('/_autocomplete', methods=['GET'])
def autocomplete():
    return Response(json.dumps(genes), mimetype='application/json')

@app.route('/_exportGeneList', methods=['POST'])
def export():
    form = SearchForm(request.form)
    fileName = 'exportedGeneList.txt'
    if request.method == 'POST':
        f = open(fileName, 'w')
        f.write('%s' % genelist[:-1])
        f.close()
        print(genelist)
    return ('', 204)
#    return render_template('index.html', form=form)

@app.route('/_downloadData', methods=['POST'])
def createUrl():
    form = SearchForm(request.form)
    if request.method == 'POST':
        print(genelist)
        url = "http://api.brain-map.org/api/v2/data/query.json?criteria=service::human_microarray_expression[probes$in"
        url += genelist
        url = url[:-1]
        url += "][donors$eq"
        donorIds = ['15496','14380','15697','9861','12876','10021']
        for d in donorIds:
            readapi.queryAPI(url, d)
        return ('', 204)
        #return render_template('index.html', form=form)

if __name__ == '__main__':
    app.run(debug=True)
