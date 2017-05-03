from flask import Flask, render_template, request, jsonify, json
import backend
# Initialize the Flask application
app = Flask(__name__)

@app.route('/')
def index():
    return render_template('index.html')

@app.route('/_jugex')
def jugex():
    res = backend.performJugex()
    return jsonify(result=res)

if __name__ == '__main__':
    debug = true
