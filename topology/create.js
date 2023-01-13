const xlsx = require('node-xlsx').default;
const classifyPoint = require("robust-point-in-polygon");
const axios = require('axios').default;
const fs = require('fs');




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
        return {valid: true, duration: minutes};
    }
    return {valid: false, duration:0};
}



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
    let pointId = 1;
    for(const [k, v] of Object.entries(points)) {
        for(const point of v){
            const timeWindow = TIME_WINDOWS[getRandom(0, 8)];
            const visitTime = VISIT_TIMES[getRandom(0, 8)];
            const profit = getRandom(5, 20);
            const lat = (+point.lat).toFixed(5);
            const lon = (+point.lon).toFixed(5);
            const sight = `${pointId++} ${lat} ${lon} ${profit} ${visitTime} ${timeWindow.startTime} ${timeWindow.endTime}`;
            console.log(`Writing node: ${sight}`);
            writeStream.write(`${sight}\n`);
        }

    }
}

const writeRoutes = async(points, writeStream) => {
    let pool = [];
    for(const [k, v] of Object.entries(points)){
        pool.push(...v);
    }
    let routeId = 0;
    for(const [i, pointX] of pool.entries()){
        for(const [j,pointY] of pool.entries()){
            const {valid, duration} = await getRoute(pointX, pointY);
            if(!valid){
                console.log(`invalid route for points ${JSON.stringify(pointX)} and ${JSON.stringify(pointY)}`);
                process.exit()
            }
            const route = `D ${routeId++} ${i} ${j} ${duration}`;
            writeStream.write(`${route}\n`)
        }
    }

    // const file = fs.createWriteStream('AthensRoutes.txt');
    // file.on('error', (error) => { console.log(error) });
    // routes.forEach(v => { file.write(v.join(', ') + '\n'); });
    // file.end();
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

    const writeStream = fs.createWriteStream('AthensTopology.txt');
    const pathName = writeStream.path;
    const data = fs.readFileSync('./athens_box.geojson', 'utf8');
    const geojson = JSON.parse(data);
    const polygon = geojson.features[0].geometry.coordinates[0];
    let pois = {};
    const validPois = {}; 
    let totalValidPoints = 0;

    const dir = "./pois";
    fs.readdirSync(dir).forEach(filename => {
        const content = xlsx.parse(`${dir}/${filename}`);
        let points = content[0].data.map(x => {return {lat: x[3], lon: x[2]}});
        points = points.slice(0, 2);
        pois[filename.split(".")[0]] = points;
    });

    for(const [k, v] of Object.entries(pois)) {
        validPois[k] = v.filter(x => {
            return inside(x, polygon);
        });
        totalValidPoints += validPois[k].length;
    }
    
    
    writeStream.on('finish', () => {
        console.log(`wrote all routes to file ${pathName}`);
    });
    writeStream.on('error', (err) => {
        console.error(`Error while writing routes: ${pathName} => ${err}`)
    });

    writeStream.write(`${totalValidPoints} ${totalValidPoints * totalValidPoints}\n`);

    await writeSights(validPois, writeStream);

    await writeRoutes(validPois, writeStream);
    
    writeStream.end();
}


create();






