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

const readFile = async() => {
    
}

const getRoute = async({lat:latX, lon:lonX}, {lat:latY, lon:lonY}) => {
    const res = await axios.get('http://localhost:8080/otp/routers/default/plan',{
        params: {
            fromPlace:`${lonX},${latX}`,
            toPlace:`${lonY},${latY}`,
            time:'7:54pm',
            date:'12-23-2022',
            mode:'CAR_PICKUP',
            arriveBy:'false',
            wheelchair:'false',
            debugItineraryFilter:'false',
            locale:'en'
        }
    })
    const {data:{plan: {itineraries}}} = res;
    if(itineraries?.length > 0){
        console.log(res);
        const duration = itineraries[0]?.duration;
        return duration;
    }
    return 0;
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

const writeSights = (validPoints, writeStream) => {
    for(const [i,point] of validPoints.entries()){
        const timeWindow = TIME_WINDOWS[getRandom(0, 8)];
        const visitTime = VISIT_TIMES[getRandom(0, 8)];
        const profit = getRandom(5, 20);
        const lat = (+point.lat).toFixed(5);
        const lon = (+point.lon).toFixed(5);
        const sight = `${i}\t${lat}\t${lon}\t${profit}\t${visitTime}\t${timeWindow.startTime}\t${timeWindow.endTime}`;
        console.log(`Writing node: ${sight}`);
        writeStream.write(`${sight}\n`);
    }
}

const writeRoutes = async(validPoints, writeStream) => {
    
    for(const [i, pointX] of validPoints.entries()){
        for(const [j,pointY] of validPoints.entries()){
            const duration = await getRoute(pointX, pointY);
            const route = `${i+j}\t${i}\t${j}\t${duration}\t0\t0\t1439`;
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
    let validPois = {};
    let totalValidPoints = 0;

    const files = [];
    const dir = "./pois";
    fs.readdirSync(dir).forEach(filename => {
        const content = xlsx.parse(`${dir}/${filename}`);
        let points = content[0].data.map(x => {return {lon: x[2], lat: x[3]}});
        points = points.slice(0, 5);
        console.log(`pois length: ${points.length}`);
        pois[filename.split(".")[0]] = points;
    
        // if (isFile) files.push({ filepath, name, ext, stat });
    });
    

    for(const [k, v] of Object.entries(pois)) {
        const validPoints = v.filter(x => {
            return inside(x, polygon);
        });
        validPois[k] = validPoints;
        totalValidPoints += validPoints.length;
    }
    console.log(totalValidPoints);
    
    writeStream.on('finish', () => {
        console.log(`wrote all routes to file ${pathName}`);
    });
    writeStream.on('error', (err) => {
        console.error(`Error while writing routes: ${pathName} => ${err}`)
    });



    writeStream.write(`${totalValidPoints} ${totalValidPoints * totalValidPoints}\n`)

    for(const [k, v] of Object.entries(pois)) {
        await writeSights(v, writeStream);
    }

    for(const [k, v] of Object.entries(pois)) {
        await writeRoutes(v, writeStream);
    }
    
    writeStream.end();
}


create();






