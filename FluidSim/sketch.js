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
/*
class Grid {
    constructor() {
      this.ylen = 12;
      this.xlen = 12;
      this.cells = ylen * xlen;
      this.u = new Float32Array(this.cells); // horizontal velo at each location
      this.v = new Float32Array(this.cells); // vertical velo at each location
      this.newU = new Float32Array(this.cells); //  temp storage for horizaontal velo
      this.newV = new Float32Array(this.cells); // temp storage for vertical velo
      this.p = new Float32Array(this.cells); // pressure at each location
      this.s = new Float32Array(this.cells); // solid or not at each location (0 is solid)
      this.m = new Float32Array(this.cells); // mass at each location
      this.newM = new Float32Array(this.cells); // temp storage for masses
      this.m.fill(1.0);
    }
    
    integrate(dt, grav) {
        //assumes 1 cell buffer on outside of grid
        for (let i = 1; i <= this.xlen-1; i++) {
            for (let j = 1; j < this.ylen-1; j++) {
                if (this.s[i * ylen + j] != 0.0 && this.s[i * ylen + j - 1] != 0.0) {
                    //^^ checks if this node is on the ground (s is 1d list) (v is 1d list)
                    let v_old = this.v[i * ylen + j];
                    this.v[i * ylen + j] = v_old + grav*dt; // basic euler estimation
                }
            }
        }
    }
    
  }*/