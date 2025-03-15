let grid = [];
let size = 50;
let u = [];
let v = [];
function setup() 
{
    createCanvas(800, 600);
    background(0);
    grid = Array.from({ length: displayWidth/size }, () => Array(displayHeight/size).fill(0))
    u = grid.length;
    v = grid[0].length;
    initCells();
    drawCells(grid);
}

function draw()
{
    drawCells(grid);
}

function mouseDragged() {
    let col = Math.floor(mouseX/ size);
    let row = Math.floor(mouseY/ size);
    grid[col][row].value = 255; // Set the cell value to white
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

class Cell{
    constructor(x, y, size, value){
        this.x = x;
        this.y = y;
        this.size = size;
        this.value = value;
    }
}

