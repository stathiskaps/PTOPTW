const xlsx = require('node-xlsx').default;
const classifyPoint = require("robust-point-in-polygon");
const axios = require('axios').default;
const fs = require('fs');

var topology = {
    nodes: [],
    routes: [],
};

class BoundingBox {
    constructor(botLat, topLat, leftLon, rightLon){
        this.botLat = botLat;
        this.topLat = topLat;
        this.leftLon = leftLon;
        this.rightLon = rightLon;
    }

    isInside({lat, lon}){
        if (lon >= this.leftLon && lon <= this.rightLon && lat >= this.botLat && lat <= this.topLat){
            return true;
        }
        return false;
    }
}

getRandom = (min, max) => { 
    return Math.floor(Math.random() * (max - min + 1) + min);
}

var found = 0;
var notFound = 0;

const getRoute = async({lat:latX, lon:lonX}, {lat:latY, lon:lonY}) => {

    if(latX == latY && lonX == lonY) {
        return {valid: true, duration: 0};
    }

    const params = {
        fromPlace:`${latX}, ${lonX}`,
        toPlace:`${latY}, ${lonY}`,
        time:'7:54pm',
        date:'12-31-2022',
        mode:'CAR_PICKUP',
        arriveBy:'false',
        wheelchair:'false',
        debugItineraryFilter:'false',
        locale:'en'
    };

    console.log(params);

    const res = await axios.get('http://localhost:8080/otp/routers/default/plan',{ params });

    // console.log(res)
    const {data:{plan: {itineraries}}} = res;
    if(itineraries?.length > 0){
        const minutes = Math.ceil(itineraries[0]?.duration/60);
        const steps = itineraries[0]?.legs[0]?.steps;
        const path = steps.map(s => [s.lat, s.lon]);
        return {valid: true, duration: minutes, path};
    }
    return {valid: false, duration:0, path: []};
}

const TimeBudget = {openTime: 0, closeTime: 760};

const TIME_WINDOWS = [
    {startTime: 540, endTime: 960},
    {startTime: 720, endTime: 900},
    {startTime: 600, endTime: 1200},
    {startTime: 720, endTime: 1080},
    {startTime: 600, endTime: 1080},
    {startTime: 400, endTime: 700},
    {startTime: 300, endTime: 770},
    {startTime: 200, endTime: 700},
    {startTime: 550, endTime: 840},
];

const VISIT_TIMES = [20, 23, 14, 13, 17, 30, 19, 10, 18];

const SCORES = [6, 7, 10, 12, 13, 16, 18, 19, 21];

const writeSights = (points, writeStream) => {
    let pointId = 0;
   
    const depotPoint = points[0];
    const poiPoints = points.slice(1);

    const depot = `${pointId++} ${depotPoint.lat} ${depotPoint.lon} 0 0 0 0 ${TimeBudget.openTime} ${TimeBudget.closeTime}`;
    console.log(`Writing depot: ${depot}`);
    writeStream.write(`${depot}\n`); 

    for(const point of poiPoints) {
        const timeWindow = TIME_WINDOWS[getRandom(0, 8)];
        const visitTime = VISIT_TIMES[getRandom(0, 8)];
        const profit = SCORES[getRandom(0, 8)];
        const lat = (+point.lat).toFixed(5);
        const lon = (+point.lon).toFixed(5);
        const sight = `${pointId++} ${lat} ${lon} ${visitTime} ${profit} 0 0 ${timeWindow.startTime} ${timeWindow.endTime}`;

        console.log(`Writing node: ${sight}`);
        writeStream.write(`${sight}\n`);

        const node = {id: pointId, lat: parseFloat(lat), lon: parseFloat(lon), visit_time: visitTime, profit, time_window:{start_time: timeWindow.startTime, end_time: timeWindow.endTime}};
        topology.nodes.push(node);
    }
}

const writeRoutes = async(points, writeStream) => {
    let routeId = 0;
    for(let i = 0; i < points.length; ++i){
        for(let j = 0; j < points.length; j++){
            const {valid, duration, path} = await getRoute(points[i], points[j]);
            if(!valid){
                console.log(`invalid route for points ${JSON.stringify(points[i])} and ${JSON.stringify(points[j])}`);
                process.exit()
            }
            const route = `D ${routeId++} ${i} ${j} ${duration}`;

            const r = {id: routeId, from: i, to: j, duration, path};
            topology.routes.push(r);
            writeStream.write(`${route}\n`);
        }
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

const create = async () => {
    const writeStream = fs.createWriteStream('../build/instances/Custom/AthensTopology.txt');
    const pathName = writeStream.path;
    const data = fs.readFileSync('./athens_box.geojson', 'utf8');
    const geojson = JSON.parse(data);
    const polygon = geojson.features[0].geometry.coordinates[0];
    let totalPoints = {};
    let totalValidPoints = [{lat: 37.97616, lon: 23.73530}]; 

    const dir = "./pois";
    fs.readdirSync(dir).forEach(filename => {
        const content = xlsx.parse(`${dir}/${filename}`);
        let points = content[0].data.map(x => {return {lat: x[3], lon: x[2]}});
        points = points.slice(0, 1);
        totalPoints[filename.split(".")[0]] = points;
    });

    for(const [k, v] of Object.entries(totalPoints)) {
        const validPoints = v.filter(x => {
            return inside(x, polygon);
        });
        totalValidPoints.push(...validPoints);
    }
    
    writeStream.on('finish', () => {
        console.log(`wrote all routes to file ${pathName}`);
    });
    writeStream.on('error', (err) => {
        console.error(`Error while writing routes: ${pathName} => ${err}`)
    });

    writeStream.write(`${totalValidPoints.length} ${totalValidPoints.length * totalValidPoints.length}\n`);
    writeStream.write(`0 0\n`);

    await writeSights(totalValidPoints, writeStream);
    
    await writeRoutes(totalValidPoints, writeStream);

    console.log(topology);

    const outputFilename = "./topology.json";
    fs.writeFile(outputFilename, JSON.stringify(topology, null, 4), function(err) {
        if(err) {
          console.log(err);
        } else {
          console.log("JSON saved to " + outputFilename);
        }
    }); 
    
    writeStream.end();
}


create();






