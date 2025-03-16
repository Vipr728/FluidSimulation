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

function UpdateCells(u, v, grid){
    let newGrid= [];
    for(let x = 0; x< u; x++){ 
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