
const colors = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f"];

const categories = {
    Hotel: {prefix: "fas", icon:"fa-hotel"},
    ArtLife: {prefix: "fas", icon: "fa-palette"},
    SightSeeings: {prefix: "fas", icon:"fa-monument"},
    // Casino: {prefix: "fas", icon: "fa-coins"},
    GreekFood: {prefix: "fas", icon: "fa-utensils"},
    ShoppingCenters: {prefix: "fas", icon: "fa-store"},
    Theaters: {prefix: "fas", icon: "fa-theater-masks"},
    VillageCinemas: {prefix: "fas", icon: "fa-video"},
    Museums: {prefix: "fas", icon: "fa-landmark"},
};

function compare( a, b ) {
    if ( a.profit < b.profit ){
      return 1;
    }
    if ( a.profit > b.profit ){
      return -1;
    }
    return 0;
}

const rowsHTML = () => {
    const categoriesArray = Object.entries(categories).map(([name, { prefix, icon, profitRange, count }]) => ({ name, prefix, icon, profitRange, count }));

    console.log(categoriesArray)
    categoriesArray.sort(compare);
    return categoriesArray.map(v =>`
        <tr><td><i class='${v.prefix} ${v.icon}'></i></td><td>${v.name}</td><td>[${v.profitRange || '0'}]</td><td>${v.count}</td></tr>`).join('')
};

// Create a map instance
var map = L.map('mapid').setView([37.9786, 23.7274], 13);

// Add a tile layer to the map
L.tileLayer('https://{s}.tile.openstreetmap.org/{z}/{x}/{y}.png', {
    attribution: 'Map data &copy; <a href="https://www.openstreetmap.org/">OpenStreetMap</a> contributors',
    maxZoom: 18,
}).addTo(map);

const drawVisited = (ta, color, categoryName) => {
    let customMarker = L.ExtraMarkers.icon({
        icon: categories[categoryName]?.icon,
        markerColor: color,
        svg: true,
        shape: 'square',
        prefix: categories[categoryName]?.prefix
    });

    const point = [ta.lat, ta.lon];
    const html = `
        <table class="table">
        <tbody>
            <tr><th scope="row">Profit</th><td>${ta.profit}</td></tr>
            <tr><th scope="row">Arrival Time</th><td>${ta.arrival_time}</td></tr>
            <tr><th scope="row">Wait duration</th><td>${ta.wait_duration}</td></tr>
            <tr><th scope="row">Start of visit time</th><td>${ta.start_visit_time}</td></tr>
            <tr><th scope="row">Departure time</th><td>${ta.departure_time}</td></tr>
            <tr><th scope="row">Shift</th><td>${ta.shift}</td></tr>
            <tr><th scope="row">Time Window</th><td>[${ta.time_window.open}-${ta.time_window.close}]</td></tr>
        </tbody>
        </table>
    `;
    const marker = L.marker(point, { icon: customMarker }).addTo(map);
    marker.on('click', function(e) {
        // do something when the marker is clicked
        const popup = L.popup()
            .setLatLng(e.latlng)
            .setContent(html)
            .openOn(map);
    });

}

const readFile = async () => {
    let response = await fetch('http://localhost:8001/solution.json');
    const solution = await response.json();

    response = await fetch('http://localhost:8001/topology.json');
    const topology = await response.json();

    

    console.log("topology", topology);

    for(const [k, v] of Object.entries(topology.categories)){
        if(categories[k]){
            categories[k]["profitRange"] = v["profit_range"];
            categories[k]["count"] = v["count"];
        }
    }

    console.log(categories)

    if(solution.unvisited){
        for(const ta of solution.unvisited){
            const node = topology.nodes.find(x => x.id == ta.id);
            const customMarker = L.ExtraMarkers.icon({
                icon: categories[node.category]?.icon,
                markerColor: "#8E8E8E",
                svg: true,
                shape: 'square',
                prefix: categories[node.category]?.prefix
            });
            
            const point = [ta.lat, ta.lon];
            const html = `
                <table class="table">
                    <tbody>
                        <tr><th scope="row">Profit</th><td>${ta.profit}</td></tr>
                        <tr><th scope="row">Visit duration</th><td>${ta.visit_duration}</td></tr>
                        <tr><th scope="row">Time Window</th><td>[${ta.time_window.open}-${ta.time_window.close}]</td></tr>
                    </tbody>
                </table>
            `;
            const marker = L.marker(point, { icon: customMarker} ).addTo(map);
            marker.on('click', function(e) {
                const popup = L.popup()
                    .setLatLng(e.latlng)
                    .setContent(html)
                    .openOn(map);
            });
        }
    }

    if(solution.walks){
        for(const [i, walk] of solution.walks.entries()){
            for(const ta of walk){
                let node;
                if(ta.id == 0){
                    // depots.push(ta);
                    drawVisited(ta, 'black', 'Hotel');
                } else {
                    node = topology.nodes.find(x => x.id == ta.id);
                    drawVisited(ta, colors[i], node.category);
                }
            }
        }
        
    }

    var polylineGroup = L.featureGroup();
    if(solution.walks){
        for(const [i, walk] of solution.walks.entries()){
            if(walk.length <= 2) {
                continue;
            }
            let pathPoints = [];
            for(let j = 0; j < walk.length-1; j+=1){
                const fromPoint = walk[j];
                const toPoint = walk[j+1];
                
                const fromId = fromPoint.id;
                const toId = toPoint.id;

                const path = topology.routes.find(route => route.from == fromId && route.to == toId)?.path;
                if(path){
                    pathPoints.push(...path);
                }
            }
            var polyline = L.polyline(pathPoints, {color: colors[i]});
            polylineGroup.addLayer(polyline);
            // Zoom the map to fit the polyline
            // map.fitBounds(polyline.getBounds());
        }
       
    }

    polylineGroup.addTo(map);
    map.fitBounds(polylineGroup.getBounds());

    var legend = L.control({ position: 'topleft' });

    // add the onAdd method to the control
    legend.onAdd = function(map) {
        // create the container element
        var container = L.DomUtil.create('div', 'leaflet-control legend');

        // add the text to the container
        container.innerHTML = `
            <table class="table">
                <thead>
                    <tr><th>Icon</th><th>Meaning</th><th>Profit range</th><th>Count</th></tr>
                </thead>
                <tbody>
                ${rowsHTML()}
                </tbody>
            </table>
            <div>Score=${solution.score}</div>
        `;

        // return the container element
        return container;
    };

    legend.addTo(map);

}

readFile();
