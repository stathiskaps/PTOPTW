
const colors = ["#2A81CB", "#9C2BCB", "#CB2B3E", "#2AAD27"];
const icons = [
    "marker-icon-blue.png",
    "marker-icon-violet.png",
    "marker-icon-red.png",
    "marker-icon-green.png"
];

const extraColors = ["purple", "red", "orange", "green"];

const categories = {
    ArtLife: {prefix: "fas", icon: "fa-palette"},
    Casino: {prefix: "fas", icon: "fa-coins"},
    GreekFood: {prefix: "fas", icon: "fa-utensils"},
    ShoppingCenters: {prefix: "fas", icon: "fa-store"},
    SPA: {prefix: "fas", icon: "fa-hot-tub"},
    Theaters: {prefix: "fas", icon: "fa-theater-masks"},
    VillageCinemas: {prefix: "fas", icon: "fa-video"},
    Museums: {prefix: "fas", icon: "fa-landmark"},
    LunaPark: {prefix: "fas", icon: "fa-tree"},
};

// Create a map instance
var map = L.map('mapid').setView([37.9786, 23.7274], 13);

// Add a tile layer to the map
L.tileLayer('https://{s}.tile.openstreetmap.org/{z}/{x}/{y}.png', {
    attribution: 'Map data &copy; <a href="https://www.openstreetmap.org/">OpenStreetMap</a> contributors',
    maxZoom: 18,
}).addTo(map);

const readFile = async () => {
    let response = await fetch('http://localhost:8001/solution.json');
    const solution = await response.json();

    response = await fetch('http://localhost:8001/topology.json');
    const topology = await response.json();

    console.log("topology", topology);

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
            <div class="metadata">
            <div>p: ${ta.profit}</div>
            <div>tw: [${ta.time_window.open} - ${ta.time_window.close}]</div>
            </div>
            `;
            L.marker(point, { icon: customMarker} ).addTo(map);
            // L.marker(point, { icon: L.divIcon({
            //     html: html,
            //     className: 'text-below-marker',
            // })}).addTo(map);
        }
    }

    if(solution.walks){
        for(const [i, walk] of solution.walks.entries()){
            
            let icon = new L.Icon({
                iconUrl: `https://raw.githubusercontent.com/pointhi/leaflet-color-markers/master/img/${icons[i]}`,
                shadowUrl: 'https://cdnjs.cloudflare.com/ajax/libs/leaflet/0.7.7/images/marker-shadow.png',
                iconSize: [25, 41],
                iconAnchor: [12, 41],
                popupAnchor: [1, -34],
                shadowSize: [41, 41]
            });

           console.log(topology.nodes);

            for(const ta of walk){

                let customMarker, node;
                if(ta.id == 0){
                    customMarker = L.ExtraMarkers.icon({
                        icon: "fa-star",
                        markerColor: colors[i],
                        svg: true,
                        shape: 'square',
                        prefix: "fas"
                    });
                } else {
                    node = topology.nodes.find(x => x.id == ta.id);
                    customMarker = L.ExtraMarkers.icon({
                        icon: categories[node.category]?.icon,
                        markerColor: colors[i],
                        svg: true,
                        shape: 'square',
                        prefix: categories[node.category]?.prefix
                    });
                }  
                
                const point = [ta.lat, ta.lon];
                const html = `
                <div style="background:${colors[i]}" class="metadata in">
                <div>${node?.category}</div>
                <div>p: ${ta.profit}</div>
                <div>tw: [${ta.time_window.open} - ${ta.time_window.close}]</div>
                </div>
                `;
                L.marker(point, { icon: customMarker }).addTo(map);
                // L.marker(point, { icon: L.divIcon({
                //     html: html,
                //     className: 'text-below-marker',
                // })}).addTo(map);
            }
        }
    }

    const hotel = solution.walks[0][0];
    console.log("hotel", hotel)

    console.log("solution", solution)

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
            var polyline = L.polyline(pathPoints, {color: colors[i]}).addTo(map);
            // Zoom the map to fit the polyline
            // map.fitBounds(polyline.getBounds());
        }
       
    }

    // create the memorandum control
    var legend = L.control({ position: 'topleft' });

    // add the onAdd method to the control
    legend.onAdd = function(map) {
        // create the container element
        var container = L.DomUtil.create('div', 'leaflet-control legend');

        // add the text to the container
        container.innerHTML = `
            <table class="table">
                <thead>
                    <tr><th>Icon</th><th>Meaning</th><th>Score</th></tr>
                </thead>
                <tbody>
                    <tr><td><i class='fas fa-palette'></i></td><td>ArtLife</td><td>30</td></tr>
                    <tr><td><i class='fas fa-coins'></i></td><td>Casino</td><td>30</td></tr>
                    <tr><td><i class='fas fa-utensils'></i></td><td>Greek Food</td><td>30</td></tr>
                    <tr><td><i class='fas fa-store'></i></td><td>Shopping Center</td><td>30</td></tr>
                    <tr><td><i class='fas fa-hot-tub'></i></td><td>Spa</td><td>30</td></tr>
                    <tr><td><i class='fas fa-theater-masks'></i></td><td>Theater</td><td>30</td></tr>
                    <tr><td><i class='fas fa-video'></i></td><td>Cinema</td><td>30</td></tr>
                    <tr><td><i class='fas fa-landmark'></i></td><td>Museum</td><td>30</td></tr>
                    <tr><td><i class='fas fa-tree'></i></td><td>Luna Park</td><td>30</td></tr>
                </tbody>
            </table>
        `;
        // container.innerHTML = "<ul>" +
        // "<li><i class='fas fa-palette'></i><span>ArtLife</span></li>" +
        // "<li><i class='fas fa-coins'></i><span>Casino</span></li>" +
        // "<li><i class='fas fa-utensils'></i><span>Greek Food</span></li>" +
        // "<li><i class='fas fa-store'></i><span>Shopping Center</span></li>" +
        // "<li><i class='fas fa-hot-tub'></i><span>Spa</span></li>" +
        // "<li><i class='fas fa-theater-masks'></i><span>Theater</span></li>" +
        // "<li><i class='fas fa-video'></i><span>Cinema</span></li>" +
        // "<li><i class='fas fa-landmark'></i><span>Museum</span></li>" +
        // "<li><i class='fas fa-tree'></i><span>Luna Park</span></li>" +
        // "</ul>";

        // return the container element
        return container;
    };

    legend.addTo(map);

}

readFile();
