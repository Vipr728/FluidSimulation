let grid = [];
let size = 20;
let u = [];
let v = [];
let add = false;
let svgWidth = 20;
let svgHeight = 20;
let svg = [];
function setup() 
{
    createCanvas(800, 600);
    background(0);
    grid = Array.from({ length: width/size }, () => Array(height/size).fill(0))
    u = grid.length;
    v = grid[0].length;
    svg = loadStrings("arrow.txt");
    initCells();
    drawCells(grid);
    let toggle = createButton("Add Circle")
    toggle.mousePressed(function() {
        add = !add;});
}

function draw()
{
    drawCells(grid);
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

function drawCells(grid){
    for(row of grid){
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

