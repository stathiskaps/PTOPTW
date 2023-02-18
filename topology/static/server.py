import http.server
import socketserver

class MyRequestHandler(http.server.SimpleHTTPRequestHandler):
    def end_headers(self):
        self.send_header('Access-Control-Allow-Origin', '*')
        super().end_headers()

with socketserver.TCPServer(("", 8001), MyRequestHandler) as httpd:
    print("serving at port 8001")
    httpd.serve_forever()
