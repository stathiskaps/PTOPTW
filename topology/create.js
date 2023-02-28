const xlsx = require('node-xlsx').default;
const classifyPoint = require("robust-point-in-polygon");
const polyline = require('@mapbox/polyline');
const axios = require('axios').default;
const fs = require('fs');

var topology = {
    categories: {Hotel: {profit: 0, count: 1}},
    nodes: [],
    routes: [],
};

const TimeBudget = {openTime: 0, closeTime: 720};

const TIME_WINDOWS = [
    {startTime: 0, endTime: 180},
    {startTime: 0, endTime: 360},
    {startTime: 60, endTime: 300},
    {startTime: 180, endTime: 540},
    {startTime: 360, endTime: 720},
    {startTime: 540, endTime: 720},
    {startTime: 540, endTime: 900},
];

const VISIT_TIME_RANGES = [[10, 15], [20, 25], [30, 35], [40, 45]];

const PROFIT_RANGES = [[5, 10], [10, 15], [15, 20], [20, 25], [25, 30], [30, 40], [40, 50], [50, 60], [60, 70]];

getRandom = ([min, max]) => { 
    return Math.floor(Math.random() * (max - min + 1) + min);
}

const getRoute = async({lat:latX, lon:lonX}, {lat:latY, lon:lonY}) => {

    if(latX == latY && lonX == lonY) {
        return {valid: true, duration: 0};
    }

    const params = {
        fromPlace:`${latX}, ${lonX}`,
        toPlace:`${latY}, ${lonY}`,
        time:'7:54pm',
        date:'12-31-2023',
        mode:'CAR_PICKUP',
        arriveBy:'false',
        wheelchair:'false',
        debugItineraryFilter:'false',
        locale:'en'
    };

    const res = await axios.get('http://localhost:8080/otp/routers/default/plan',{ params });
    const {data:{plan: {itineraries}}} = res;
    if(itineraries?.length > 0){
        const leg = itineraries[0].legs[0]; // assume we want the first leg of the itinerary
        const encodedPolyline = leg.legGeometry.points; // get the encoded polyline from the leg
        const path = polyline.decode(encodedPolyline); // decode the polyline to get the coordinates
        const minutes = Math.ceil(itineraries[0]?.duration/60);
        return {valid: true, duration: minutes, path};
    }
    return {valid: false, duration:0, path: []};
}

const addRoutes = async(points) => {
    let routeId = 0;
    let counter = 0, total = points.length*points.length;
    for(let i = 0; i < points.length; ++i){
        for(let j = 0; j < points.length; j++){
            console.log(`Requesting route from ${JSON.stringify(points[i])} to ${JSON.stringify(points[j])}`)
            const {valid, duration, path} = await getRoute(points[i], points[j]);
            if(!valid){
                console.log(`invalid route for points ${JSON.stringify(points[i])} and ${JSON.stringify(points[j])}`);
                process.exit()
            }
            console.log(`Got route ${counter+=1} of ${total}`);
            const r = {id: routeId, from: i, to: j, duration, path};
            topology.routes.push(r);
        }
    }
}



const addSights = (points) => {
    let pointId = 0;
   
    const depotPoint = points[0];
    const pois = points.slice(1);

    const depot = {id: pointId, lat: parseFloat(depotPoint.lat), lon: parseFloat(depotPoint.lon), visit_time: 0, category: "Hotel", time_window:{start_time: TimeBudget.openTime, end_time: TimeBudget.closeTime}};
    topology.nodes.push(depot);

    for(const point of pois) {
        let visitTime, timeWindow;
        if(pointId < pois.length * 0.25){
           visitTime = getRandom(VISIT_TIME_RANGES[0]);
        } else if(pointId < pois.length * .5){
            visitTime = getRandom(VISIT_TIME_RANGES[1]);
        } else if ([pointId < pois.length * 0.75]){
            visitTime = getRandom(VISIT_TIME_RANGES[2]);
        } else {
            visitTime = getRandom(VISIT_TIME_RANGES[3]);
        }

        if(pointId < pois.length * 0.14){
            timeWindow = TIME_WINDOWS[0];
        } else if (pointId < pois.length * 0.28){
            timeWindow = TIME_WINDOWS[1];
        } else if (pointId < pois.length * 0.42){
            timeWindow = TIME_WINDOWS[2];
        } else if (pointId < pois.length * 0.56){
            timeWindow = TIME_WINDOWS[3];
        } else if (pointId < pois.length * 0.70){
            timeWindow = TIME_WINDOWS[4];
        } else if (pointId < pois.length * 0.84){
            timeWindow = TIME_WINDOWS[5];
        } else {
            timeWindow = TIME_WINDOWS[6];
        }
        const lat = (+point.lat).toFixed(5);
        const lon = (+point.lon).toFixed(5);        
        const node = {id: ++pointId, lat: parseFloat(lat), lon: parseFloat(lon), visit_time: visitTime, category: point.category, time_window:{start_time: timeWindow.startTime, end_time: timeWindow.endTime}};
        topology.nodes.push(node);
    }
}

function inside(point, polygon) {
    // ray-casting algorithm based on
    // https://wrf.ecse.rpi.edu/Research/Short_Notes/pnpoly.html/pnpoly.html
    
    var x = point.lon, y = point.lat;
    
    var inside = false;
    for (var i = 0, j = polygon.length - 1; i < polygon.length; j = i++) {
        var xi = polygon[i][0], yi = polygon[i][1];
        var xj = polygon[j][0], yj = polygon[j][1];
        
        var intersect = ((yi > y) != (yj > y))
            && (x < (xj - xi) * (y - yi) / (yj - yi) + xi);
        if (intersect) inside = !inside;
    }
    
    return inside;
};

/* Randomize array in-place using Durstenfeld shuffle algorithm */
function shuffleArray(array) {
    for (var i = array.length - 1; i > 0; i--) {
        var j = Math.floor(Math.random() * (i + 1));
        var temp = array[i];
        array[i] = array[j];
        array[j] = temp;
    }
}

const create = async () => {
    const data = fs.readFileSync('./athens_box.geojson', 'utf8');
    const geojson = JSON.parse(data);
    const polygon = geojson.features[0].geometry.coordinates[0];
    let depotPoint = {lat: 37.97616, lon: 23.73530, category: "Hotel"};
    let totalPoints = [];
    let categories = [];

    const dir = "./pois";
    fs.readdirSync(dir).forEach(filename => {
        const content = xlsx.parse(`${dir}/${filename}`);
        const category = filename.split(".")[0];
        let points = content[0].data.map(x => {return {lat: x[3], lon: x[2], category }});
        points = points.slice(0, 4);
        totalPoints.push(...points);

        const randomIndex = getRandom([0, PROFIT_RANGES.length-1]);
        topology.categories[category] = {profit: getRandom(PROFIT_RANGES[randomIndex]), count: points.length};
        PROFIT_RANGES.splice(randomIndex, 1);
    });

    const validPoints = totalPoints.filter(x => {
        return inside(x, polygon);
    });

    shuffleArray(validPoints);
    validPoints.unshift(depotPoint);

    console.log(validPoints);

    await addSights(validPoints);
    await addRoutes(validPoints);

    const outputFilename = "./topology.json";
    fs.writeFile(outputFilename, JSON.stringify(topology, null, 4), function(err) {
        if(err) {
          console.log(err);
        } else {
          console.log("JSON saved to " + outputFilename);
        }
    }); 
}


create();






