from flask import Flask, render_template, request, jsonify, json, redirect, url_for, send_from_directory
from werkzeug import secure_filename
import backend, os
# Initialize the Flask application
app = Flask(__name__)
if not os.path.exists('uploads/'):
    os.makedirs('uploads/')
app.config['UPLOAD_FOLDER'] = 'uploads/'
app.config['ALLOWED_EXTENSIONS'] = set(['csv', 'gz'])
def allowed_file(filename):
    return '.' in filename and \
           filename.rsplit('.', 1)[1] in app.config['ALLOWED_EXTENSIONS']

@app.route('/')
def index():
    return render_template('index.html')

@app.route('/_jugex')
def jugex():
    res = backend.performJugex()
    return jsonify(result=res)

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


@app.route('/upload', methods=['POST'])
def upload():
    file = request.files['file']
    if file and allowed_file(file.filename):
        filename = secure_filename(file.filename)
        path = os.path.join(app.config['UPLOAD_FOLDER'], filename)
        print(path)
        file.save(path)
        extension = filename.rsplit('.', 1)[1]
        if(extension == 'csv'):
            rows = backend.readCSVFile(path)
        else:
            imgs.append(backend.readVOI(path))
        return render_template('index.html')

@app.route('/uploads/<filename>')
def uploaded_file(filename):
    return send_from_directory(app.config['UPLOAD_FOLDER'], filename)
    #backend.readCSVFile(filename)

if __name__ == '__main__':
    debug = true
