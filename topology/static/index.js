
// Create a map instance
var map = L.map('mapid').setView([37.9786, 23.7274], 13);

// Add a tile layer to the map
L.tileLayer('https://{s}.tile.openstreetmap.org/{z}/{x}/{y}.png', {
    attribution: 'Map data &copy; <a href="https://www.openstreetmap.org/">OpenStreetMap</a> contributors',
    maxZoom: 18,
}).addTo(map);

// Iterate through your array of points and add a marker for each point
var points = [[37.97616, 23.7353], [37.97964, 23.72449], [38.00183, 23.82721], [37.98310, 23.71775]];
for (var i = 0; i < points.length; i++) {
    var point = points[i];
    L.marker(point).addTo(map);
}

var circle = L.circle([37.97616, 23.7353], {
    color: 'red',
    fillColor: '#f03',
    fillOpacity: 0.5,
    radius: 500
}).addTo(map);