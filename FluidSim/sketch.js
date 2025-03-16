let grid = [];
let size = 20;
let u = [];
let v = [];
let add = false;
let svgWidth = 20;
let svgHeight = 20;
let svg = [];
let slider;
let gridInstance;
function setup() 
{
    createCanvas(800, 600);
    background(0);
    grid = Array.from({ length: width/size }, () => Array(height/size).fill(0))
    u = grid.length;
    v = grid[0].length;
    svg = loadStrings("arrow.txt");
    gridInstance = new Grid(grid,20,20);
    initCells();
    drawCells();
    slider = createSlider(0,100,9.81,0);
    slider.size(100);
    let circle = createButton('Circle');
    circle.mousePressed(() => changeObject('circle'));
    let square = createButton('Square');
    square.mousePressed(() => changeObject('square'));
    let abhi = createButton('abhi');
    abhi.mousePressed(() => changeObject('abhi'));
    function changeObject(object) {
        background(random(255)); 
        console.log('Changed to:', object);// implement
    }
    let toggle = createButton("Add Circle")
    toggle.mousePressed(function() {
        add = !add;});
}

function draw()
{
    drawCells(gridInstance.cellList);
    UpdateCells()
}

function mouseDragged() {
    let col = Math.floor(mouseX/ size);
    let row = Math.floor(mouseY/ size);
    if (col >= grid.length || row >= grid[0].length) return [-1,-1]; // Check bounds
    grid[col][row].value = 255; // Set the cell value to white
}

function mousePressed() {
    if (add) {
        parseSVG(svg[0], mouseX, mouseY);
    }
}

function drawCells(){
    for(row of gridInstance.cellList){
        for(cell of row){
            stroke( 255, 255, 255);
            fill(0,0,255, cell.value);
            rect(cell.x, cell.y, cell.size, cell.size);
        }
    }
}

function initCells(){
    for(let x = 0; x< width/size; x++){ 
        for(let y = 0; y< height/size; y++){
            grid[x][y] = new Cell(x*size, y*size, size, 0);
        } 
    }
    
}

function UpdateCells(){
    let grav = slider.value(); // use this value to have changing gravity vals
    for(let x = 0; x < gridInstance.cellList.length; x++) { 
        for(let y = 0; y < gridInstance.cellList[x].length; y++) {
            if (floor(random(0,100)) == 2) {
                gridInstance.cellList[x][y].value += 10;
                gridInstance.cellList[x][y].value = constrain(gridInstance.cellList[x][y].value, 0, 255);
            }
        }
    }
}
function parseSVG(svgString, centerX, centerY) {
    let parser = new DOMParser();
    const doc = parser.parseFromString(svgString, "image/svg+xml");
    const svgElement = doc.querySelector('svg');
    if (!svgElement) {
        console.error("Error parsing SVG: No SVG element found");
        return;
    }
    let paths = svgElement.querySelectorAll('path');
    let points = [];
    paths.forEach(path => {
        for (let i = 0; i < 1000; i++) {
            const point = path.getPointAtLength(i * 0.001 * path.getTotalLength());
            const adjustedX = (centerX-svgWidth) + point.x * size/1;
            const adjustedY = (centerY-svgWidth) + point.y * size/1;
            stroke(255, 0, 0);
            fill(255, 0, 0);
            points.push({ x: adjustedX, y: adjustedY });
            circle(adjustedX, adjustedY, 5);
            console.log(adjustedX, adjustedY);
        }
    });
    return points;
}
class Cell{
    constructor(x, y, size, value){
        this.x = x;
        this.y = y;
        this.size = size;
        this.value = value;
    }
}

class Grid {
    constructor(arr,xlen,ylen) {
        this.cellList = arr;
        this.grav = 9.81;
        this.ylen = ylen;
        this.xlen = xlen;
        this.numCells = ylen * xlen;
        this.u = new Float32Array(this.numCells); // horizontal velo at each location
        this.v = new Float32Array(this.numCells); // vertical velo at each location
        this.newU = new Float32Array(this.numCells); //  temp storage for horizaontal velo
        this.newV = new Float32Array(this.numCells); // temp storage for vertical velo
        this.p = new Float32Array(this.numCells); // pressure at each location
        this.s = new Float32Array(this.numCells); // solid or not at each location (0 is solid)
        this.m = new Float32Array(this.numCells); // mass at each location
        this.newM = new Float32Array(this.numCells); // temp storage for masses
        this.m.fill(1.0);
    }
    
    integrate(dt) {
        //assumes 1 cell buffer on outside of grid
        for (let i = 1; i <= this.xlen-1; i++) {
            for (let j = 1; j < this.ylen-1; j++) {
                if (this.s[i * ylen + j] != 0.0 && this.s[i * ylen + j - 1] != 0.0) {
                    //^^ checks if this node is on the ground (s is 1d list) (v is 1d list)
                    let v_old = this.v[i * ylen + j];
                    this.v[i * ylen + j] = v_old + this.grav*dt; // basic euler estimation
                }
            }
        }
    }
  }

class Fluid {
    constructor(density, nx, ny, h) {
        this.density = density;
        this.nx = nx+2;
        this.ny = ny+2;
        this.ncells = this.nx * this.ny;
        this.h = h;
        this.u = new Float32Array(this.ncells);
        this.v = new Float32Array(this.ncells);
        this.newU = new Float32Array(this.ncells);
        this.newV = new Float32Array(this.ncells);
        this.p = new Float32Array(this.ncells);
        this.s = new Float32Array(this.ncells);
        this.m = new Float32Array(this.ncells);
        this.newM = new Float32Array(this.ncells);
        this.m.fill(1.0);
        this.OVER_RELAXATION = 1.9;
    }

    integrate(dt, grav) {
        let n = this.ny;
        for (let i = 1; i < this.nx; i++) {
            for (let j = 1; j < this.ny-1; j++) {
                if (this.s[i*n+j] && this.s[i*n+j-1]) {
                    this.v[i*n+j] += grav*dt;
                }
            }
        }
    }

    uncompress(numIterations, dt) {
        let n = this.ny;
        let cp = this.density * this.h/dt;
        while (numIterations --> 0) {
            for (let i = 1; i < this.nx-1; i++) {
                for (let j = 1; j < this.ny-1; j++) {
                    if (!this.s[i*n+j]) continue;
                    let sx0 = this.s[(i-1)*n+j];
                    let sx1 = this.s[(i+1)*n+j];
                    let sy0 = this.s[i*n+j-1]
                    let sy1 = this.s[i*n+j+1]
                    let s = sx0 + sx1 + sy0 + sy1;
                    if (!s) continue;
                    let div = this.u[(i+1)*n+j] - this.u[i*n+j] + this.v[i*n+j+1] - this.v[i*n+j];
                    let p = -div/s * this.OVER_RELAXATION;
                    this.p[i*n+j] += cp*p;
                    
                    this.u[i*n+j] -= sx0*p;
                    this.u[(i+1)*n+j] += sx1*p;
                    this.v[i*n+j] -= sy0*p;
                    this.v[i*n+j+1] += sy1*p;
                }
            }
        }
    }

    extrapolate() {
        let n = this.ny;
        for (let i = 0; i < this.nx; i++) {
            this.u[i*n+0] = this.u[i*n+1];
            this.u[i*n+n-1] = this.u[i*n+n-2];
        }
        for (let j = 0; j < this.ny; j++) {
            this.v[0*n+j] = this.v[1*n+j];
            this.v[(this.nx-1)*n + j] = this.v[(this.nx-2)*n + j];
        }
    }

    advectVel(dt) {
        this.newU.set(this.u);
        this.newV.set(this.v);

        let n = this.ny;
        let h = this.h;
        let h2 = .5*h;

        for (let i = 1; i < this.nx; i++) {
            for (let j = 1; j < this.ny; j++) {
                // cnt++;

                if (this.s[i*n+j] && this.s[(i-1)*n+j] && j+1 < this.ny) {
                    let x = i*h;
                    let y = j*h + h2;
                    let u = this.u[i*n+j]
                    let v = this.avgV(i, j);
                    x -= dt*u;
                    y -= dt*v;
                    u = this.sampleField(x, y, 'U_FIELD');
                    this.newU[i*n+j] = u;
                }

                if (this.s[i*n+j] && this.s[i*n+j-1] && i+1 < this.nx) {
                    let x = i*h+h2;
                    let y = j*h;
                    let u = this.avgU(i, j);
                    let v = this.v[i*n+j];
                    x -= dt*u;
                    y -= dt*v;
                    v = this.sampleField(x, y, 'V_FIELD');
                    this.newV[i*n+j] = v;
                }
            }
        }

        this.u.set(this.newU);
        this.v.set(this.newV);

    }

    sampleField(x, y, field) {
        let n = this.ny;
        let h = this.h;
        let h1 = 1./h;
        let h2 = .5*h;

        x = Math.max(Math.min(x, this.nx*h), h)
        y = Math.max(Math.min(y, this.ny*h), h)

        let dx = 0.;
        let dy = 0.;

        let f ;

        switch(field) {
            case 'U_FIELD': { f = this.u; dy = h2; break}
            case 'V_FIELD': { f = this.v; dx = h2; break}
            case 'S_FIELD': {f = this.m; dx = h2; dy = h2; break}
        }

        let x0 = Math.min(Math.floor((x-dx)*h1), this.nx-1);
        let tx = ((x-dx)-x0*h)*h1;
        let x1 = Math.min(x0+1, this.nx-1)

        let y0 = Math.min(Math.floor((y-dy)*h1), this.ny-1)
        let ty = ((y-dy)-y0*h)*h1;
        let y1 = Math.min(y0+1, this.ny-1)

        let sx = 1.-tx;
        let sy = 1.-ty;

        let val = (
            sx*sy * f[x0*n+y0] + 
            tx*sy * f[x1*n+y0] +
            tx*ty * f[x1*n+y1] +
            sx*ty * f[x0*n+y1]
        );
        return val;
    }

    avgU(i, j) {
        let n = this.ny;
        let u = (
            this.u[i*n+j-1] + this.u[i*n+j] + 
            this.u[(i+1)*n+j-1] + this.u[(i+1)*n+j]
        ) * 0.25;
        return u;
    }
    
    avgV(i, j) {
        let n = this.ny;
        let v = (
            this.v[(i-1)*n+j] + this.v[i*n+j] + 
            this.v[(i-1)*n+j+1] + this.v[i*n+j+1]
        ) * .25;
        return v;
    }

    advectSmoke(dt) {
        this.newM.set(this.m);

        let n = this.ny;
        let h = this.h;
        let h2 = .5*h;
        for (let i = 1; i < this.nx-1; i++) {
            for (let j = 1; j < this.ny-1; j++) {

                if (!this.s[i*n+j]) continue;

                let u = .5 * (this.u[i*n+j] + this.u[(i+1)*n+j])
                let v = .5 * (this.v[i*n+j] + this.v[i*n+j+1])
                let x = i*h + h2 - dt*u;
                let y = j*h + h2 - dt*v;
                this.newM[i*n+j] = this.sampleField(x,y, 'S_FIELD')
            }
        }
        this.m.set(newM);
    }

    simulate(dt, grav, iterations) {
        this.integrate(dt, grav);
        this.p.fill(0.);
        this.uncompress(iterations, dt);
        
        this.extrapolate();
        this.advectVel(dt);
        this.advectSmoke(dt);

    }


}