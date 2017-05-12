from flask import Flask, render_template, request, jsonify, json, redirect, url_for, send_from_directory, Response
from werkzeug import secure_filename
import backend, readapi, os
from wtforms import TextField, Form

# Initialize the Flask application
app = Flask(__name__)
#backend.writeGeneList()
genes = backend.readGeneList()
genelist = ""

class SearchForm(Form):
    autocomp = TextField('Insert Entrez_id', id='gene_autocomplete')

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
    if request.method == 'POST':
        files = request.files.getlist("file[]")
        for file in files:
            if file and allowed_file(file.filename):
                filename = secure_filename(file.filename)
                path = os.path.join(app.config['UPLOAD_FOLDER'], filename)
                file.save(path)
        return render_template('index.html')


@app.route('/_autocomplete', methods=['GET'])
def autocomplete():
    return Response(json.dumps(genes), mimetype='application/json')

@app.route('/_downloadData', methods=['POST'])
def createUrl():
    if request.method == 'POST':
        print(genelist)
        url = "http://api.brain-map.org/api/v2/data/query.json?criteria=service::human_microarray_expression[probes$in"
        url += genelist
        url = url[:-1]
        url += "][donors$eq"
        donorIds = ['15496','14380','15697','9861','12876','10021']
        for d in donorIds:
            readapi.queryAPI(url, d)
        return ""

if __name__ == '__main__':
    app.run(debug=True)
