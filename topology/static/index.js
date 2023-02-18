
// Create a map instance
var map = L.map('mapid').setView([37.9786, 23.7274], 13);

// Add a tile layer to the map
L.tileLayer('https://{s}.tile.openstreetmap.org/{z}/{x}/{y}.png', {
    attribution: 'Map data &copy; <a href="https://www.openstreetmap.org/">OpenStreetMap</a> contributors',
    maxZoom: 18,
}).addTo(map);

const readFile = async () => {
    const response = await fetch('http://localhost:8001/solution.json');
    const solution = await response.json();
    let coords = [];

    var unvisitedIcon = new L.Icon({
        iconUrl: 'https://raw.githubusercontent.com/pointhi/leaflet-color-markers/master/img/marker-icon-grey.png',
        shadowUrl: 'https://cdnjs.cloudflare.com/ajax/libs/leaflet/0.7.7/images/marker-shadow.png',
        iconSize: [25, 41],
        iconAnchor: [12, 41],
        popupAnchor: [1, -34],
        shadowSize: [41, 41]
    });

    if(solution.unvisited){
        for(const u of solution.unvisited){
            const point = [u.lat, u.lon];
            const html = `
            <div>${u.profit} : [${u.time_window.open} - ${u.time_window.close}]</div>
            `;
            L.marker(point, { icon: unvisitedIcon} ).addTo(map);
            L.marker(point, { icon: L.divIcon({
                html: html,
                className: 'text-below-marker',
            })}).addTo(map);
        }
    }

    if(solution.walks){
        for(const walk of solution.walks){
            for(const ta of walk){
                const point = [ta.lat, ta.lon];
                coords.push(point);

                const html = `
                <div>${ta.profit} : [${ta.time_window.open} - ${ta.time_window.close}]</div>
                `;
                L.marker(point).addTo(map);
                L.marker(point, { icon: L.divIcon({
                    html: html,
                    className: 'text-below-marker',
                })}).addTo(map);
            }
        }
    }

    var polyline = L.polyline(coords, {color: '#4592DF'}).addTo(map);

    // Zoom the map to fit the polyline
    map.fitBounds(polyline.getBounds());

}

readFile();


// var circle = L.circle([37.97616, 23.7353], {
//     color: 'red',
//     fillColor: '#f03',
//     fillOpacity: 0.5,
//     radius: 500
// }).addTo(map);