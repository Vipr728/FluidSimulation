let grid = [];
    let size = 10;
    let u = [];
    let v = [];
    let add = false;
    let svgWidth = 20;
    let svgHeight = 20;
    let svg = [];
    let slider;
    let gridInstance;
    let objects = [];
    let selectedObject = null;
    let objectType = 'circle';
    let svgObjects = [];
    let svgCache = new Map();
    let svgTemplates = {};

    function preload() {
        // Load SVGs from /svgs directory
        const svgFiles = ['abhi.svg', 'cow.svg', 'george.svg'];
        svgFiles.forEach(filename => {
            // Use loadXML instead of loadStrings
            svgTemplates[filename] = loadXML(`svgs/${filename}`);
        });
    }
    
    function addSVGObject(filename) {
        const svgXML = svgTemplates[filename];
        const svgString = new XMLSerializer().serializeToString(svgXML);
        const svgData = parseSVG(svgString, 0, 0, 0.5);
        objects.push(new DraggableObject(mouseX, mouseY, 'svg', svgData));
    }

    function setup() {
        createCanvas(2000, 920);
        background(0);
        grid = Array.from({ length: width / size }, () => Array(height / size).fill(0));
        u = grid.length;
        v = grid[0].length;
        loadStrings("arrow.txt", (result) => {
            svg = result;
        });
        initCells();
        drawCells(grid);
        gridInstance = new Fluid(grid, 1000, width / size, height / size, size);
        slider = createSlider(0, 100, 9.81, 0);
        slider.size(100);
        let circle = createButton('Circle');
        circle.mousePressed(() => changeObject('circle'));
        let square = createButton('Square');
        square.mousePressed(() => changeObject('square'));
        let abhi = createButton('abhi');
        abhi.mousePressed(() => changeObject('abhi'));
        function changeObject(object) {
            background(random(255));
            console.log('Changed to:', object); // implement
        }
        let toggle = createButton("Add Circle");
        toggle.mousePressed(function () {
            add = !add;
        });
        let spoutSlider = createSlider(0,1000, 250, 1);
        spoutSlider.input(() => {
            gridInstance.spoutIntensity = spoutSlider.value();
        });
        let addCircle = createButton('Add Circle');
        addCircle.mousePressed(() => addNewObject('circle'));
        let addSquare = createButton('Add Square');
        addSquare.mousePressed(() => addNewObject('square'));

        function addNewObject(type) {
        objects.push(new DraggableObject(width/2, height/2, 60, type));
        }
        Object.keys(svgTemplates).forEach(name => {
            let btn = createButton(`Add ${name}`);
            btn.mousePressed(() => addSVGObject(name));
          });
    }

    function draw() {
        background(0);
        var r = frameCount % 200 * Math.sqrt(2);
        ellipse(100, 100, r, r);
        // Draw fluid first
        drawCells(gridInstance.cellList);
        gridInstance.simulate(0.1, 9.81, 20);
        
        // Then draw objects on top
        for(let obj of objects) {
            obj.update();
            obj.show();
        }
    }

    function mouseDragged() {
        let col = Math.floor(mouseX / size);
        let row = Math.floor(mouseY / size);
        if (col >= grid.length || row >= grid[0].length) return; // Check bounds
        grid[col][row].value = 255; // Set the cell value to white
    }

    function mousePressed() {
        // Check object selection
        for(let i = objects.length-1; i >= 0; i--) {
          if(objects[i].contains(mouseX, mouseY)) {
            selectedObject = objects[i];
            selectedObject.dragging = true;
            return;
          }
        }
      }
    
      function mouseReleased() {
        if(selectedObject) {
          selectedObject.dragging = false;
          selectedObject = null;
        }
      }

    function drawCells(cellList) {
        for(let row of cellList) {
            for(let cell of row) {
                noStroke();
                // Use velocity-influenced color with mass-based opacity
                fill(cell.color);
                rect(cell.x, cell.y, cell.size, cell.size);
            }
        }
    }

    function initCells() {
        for (let x = 0; x < width / size; x++) {
            for (let y = 0; y < height / size; y++) {
                grid[x][y] = new Cell(x * size, y * size, size, 0);
            }
        }
    }

    // function parseSVG(svgString, centerX, centerY, scaleFactor = 1) {
    //     let points = [];
    //     let parser = new DOMParser();
    //     const doc = parser.parseFromString(svgString, "image/svg+xml");
    //     const svgElement = doc.querySelector('svg');
        
    //     // Get dimensions from viewBox if width/height not specified
    //     const viewBox = svgElement.getAttribute('viewBox')?.split(' ') || [0, 0, svgWidth, svgHeight];
    //     const svgW = parseFloat(svgElement.getAttribute('width')) || viewBox[2];
    //     const svgH = parseFloat(svgElement.getAttribute('height')) || viewBox[3];
    
    //     // Process all paths
    //     const paths = svgElement.querySelectorAll('path');
    //     paths.forEach(path => {
    //         const totalLength = path.getTotalLength();
    //         for (let i = 0; i <= 100; i++) {
    //             const point = path.getPointAtLength((i/100) * totalLength);
    //             points.push({
    //                 x: (point.x - svgW/2) * scaleFactor,
    //                 y: (point.y - svgH/2) * scaleFactor
    //             });
    //         }
    //     });
        
    //     return {
    //         points,
    //         w: svgW * scaleFactor,
    //         h: svgH * scaleFactor
    //     };
    // }
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
    /**
     * Determines if a point (x, y) is inside a given SVG object.
     *
     * @param {number} x - The x-coordinate of the point to check.
     * @param {number} y - The y-coordinate of the point to check.
     * @param {Object} obj - The SVG object to check against.
     * @param {number} obj.x - The x-coordinate of the SVG object's center.
     * @param {number} obj.y - The y-coordinate of the SVG object's center.
     * @param {Object} obj.data - The data of the SVG object.
     * @param {number} obj.data.w - The width of the SVG object.
     * @param {number} obj.data.h - The height of the SVG object.
     * @param {Array} obj.data.points - The points defining the SVG object's shape.
     * @param {number} obj.data.points[].x - The x-coordinate of a point in the SVG object's shape.
     * @param {number} obj.data.points[].y - The y-coordinate of a point in the SVG object's shape.
     * @returns {boolean} - Returns true if the point is inside the SVG object, false otherwise.
     */
    function pointInSVG(x, y, obj) {
        const localX = x - (obj.x - obj.data.w/2);
        const localY = y - (obj.y - obj.data.h/2);
        
        // Simple ray-casting check
        let inside = false;
        const points = obj.data.points;
        for (let i = 0, j = points.length - 1; i < points.length; j = i++) {
          const xi = points[i].x, yi = points[i].y;
          const xj = points[j].x, yj = points[j].y;
          
          const intersect = ((yi > localY) !== (yj > localY)) &&
            (localX < (xj - xi) * (localY - yi) / (yj - yi) + xi);
          if (intersect) inside = !inside;
        }
        return inside;
      }

    class Cell {
        constructor(x, y, size, value) {
            this.x = x;
            this.y = y;
            this.size = size;
            this.value = value;
            this.color = [0, 50, 200, 255]; // [R, G, B, A]
        }
    }

    class Fluid {
        constructor(arr, density, nx, ny, h) {
            this.cellList = arr;
            this.density = density;
            this.nx = nx + 2;
            this.ny = ny + 2;
            this.ncells = this.nx * this.ny;
            this.h = h;
            
            // Initialize all arrays
            this.u = new Float32Array(this.ncells);
            this.v = new Float32Array(this.ncells);
            this.newU = new Float32Array(this.ncells);
            this.newV = new Float32Array(this.ncells);
            this.p = new Float32Array(this.ncells);
            this.s = new Float32Array(this.ncells);
            this.m = new Float32Array(this.ncells);
            this.newM = new Float32Array(this.ncells);
    
            // Mark fluid cells (1 = fluid, 0 = solid)
            for(let i = 1; i < this.nx-1; i++) {
                for(let j = 1; j < this.ny-1; j++) {
                    this.s[i * this.ny + j] = 1;
                }
            }
            
            this.OVER_RELAXATION = 1.9;
            this.spoutIntensity = 250; // Base velocity for spout
        }
    

        integrate(dt, grav) {
            let n = this.ny;
            for (let i = 1; i < this.nx; i++) {
                // Changed loop bounds and direction
                for (let j = this.ny - 2; j >= 1; j--) {  // Start from bottom
                    if (this.s[i * n + j] && this.s[i * n + j + 1]) {  // Check cell below
                        this.v[i * n + j] += 0 * dt;
                    }
                }
            }
        }

        uncompress(numIterations, dt) {
            let n = this.ny;
            let cp = this.density * this.h / dt;
            while (numIterations-- > 0) {
                for (let i = 1; i < this.nx - 1; i++) {
                    for (let j = 1; j < this.ny - 1; j++) {
                        if (!this.s[i * n + j]) continue;
                        let sx0 = this.s[(i - 1) * n + j];
                        let sx1 = this.s[(i + 1) * n + j];
                        let sy0 = this.s[i * n + j - 1];
                        let sy1 = this.s[i * n + j + 1];
                        let s = sx0 + sx1 + sy0 + sy1;
                        if (!s) continue;
                        let div = this.u[(i + 1) * n + j] - this.u[i * n + j] + this.v[i * n + j + 1] - this.v[i * n + j];
                        let p = -div / s * this.OVER_RELAXATION;
                        this.p[i * n + j] += cp * p;

                        this.u[i * n + j] -= sx0 * p;
                        this.u[(i + 1) * n + j] += sx1 * p;
                        this.v[i * n + j] -= sy0 * p;
                        this.v[i * n + j + 1] += sy1 * p;
                    }
                }
            }
        }

        extrapolate() {
            let n = this.ny;
            for (let i = 0; i < this.nx; i++) {
                this.u[i * n + 0] = this.u[i * n + 1];
                this.u[i * n + n - 1] = this.u[i * n + n - 2];
            }
            for (let j = 0; j < this.ny; j++) {
                this.v[0 * n + j] = this.v[1 * n + j];
                this.v[(this.nx - 1) * n + j] = this.v[(this.nx - 2) * n + j];
            }
        }

        advectVel(dt) {
            this.newU.set(this.u);
            this.newV.set(this.v);

            let n = this.ny;
            let h = this.h;
            let h2 = .5 * h;

            for (let i = 1; i < this.nx; i++) {
                for (let j = 1; j < this.ny; j++) {
                    if (this.s[i * n + j] && this.s[(i - 1) * n + j] && j + 1 < this.ny) {
                        let x = i * h;
                        let y = j * h + h2;
                        let u = this.u[i * n + j];
                        let v = this.avgV(i, j);
                        x -= dt * u;
                        y -= dt * v;
                        u = this.sampleField(x, y, 'U_FIELD');
                        this.newU[i * n + j] = u;
                    }

                    if (this.s[i * n + j] && this.s[i * n + j - 1] && i + 1 < this.nx) {
                        let x = i * h + h2;
                        let y = j * h;
                        let u = this.avgU(i, j);
                        let v = this.v[i * n + j];
                        x -= dt * u;
                        y -= dt * v;
                        v = this.sampleField(x, y, 'V_FIELD');
                        this.newV[i * n + j] = v;
                    }
                }
            }

            this.u.set(this.newU);
            this.v.set(this.newV);
        }

        sampleField(x, y, field) {
            let n = this.ny;
            let h = this.h;
            let h1 = 1. / h;
            let h2 = .5 * h;

            x = Math.max(Math.min(x, this.nx * h), h);
            y = Math.max(Math.min(y, this.ny * h), h);

            let dx = 0.;
            let dy = 0.;

            let f;

            switch (field) {
                case 'U_FIELD': { f = this.u; dy = h2; break; }
                case 'V_FIELD': { f = this.v; dx = h2; break; }
                case 'S_FIELD': { f = this.m; dx = h2; dy = h2; break; }
            }

            let x0 = Math.min(Math.floor((x - dx) * h1), this.nx - 1);
            let tx = ((x - dx) - x0 * h) * h1;
            let x1 = Math.min(x0 + 1, this.nx - 1);

            let y0 = Math.min(Math.floor((y - dy) * h1), this.ny - 1);
            let ty = ((y - dy) - y0 * h) * h1;
            let y1 = Math.min(y0 + 1, this.ny - 1);

            let sx = 1. - tx;
            let sy = 1. - ty;

            let val = (
                sx * sy * f[x0 * n + y0] +
                tx * sy * f[x1 * n + y0] +
                tx * ty * f[x1 * n + y1] +
                sx * ty * f[x0 * n + y1]
            );
            return val;
        }

        avgU(i, j) {
            let n = this.ny;
            let u = (
                this.u[i * n + j - 1] + this.u[i * n + j] +
                this.u[(i + 1) * n + j - 1] + this.u[(i + 1) * n + j]
            ) * 0.25;
            return u;
        }

        avgV(i, j) {
            let n = this.ny;
            let v = (
                this.v[(i - 1) * n + j] + this.v[i * n + j] +
                this.v[(i - 1) * n + j + 1] + this.v[i * n + j + 1]
            ) * .25;
            return v;
        }

        advectSmoke(dt) {
            this.newM.set(this.m);

            let n = this.ny;
            let h = this.h;
            let h2 = .5 * h;
            for (let i = 1; i < this.nx - 1; i++) {
                for (let j = 1; j < this.ny - 1; j++) {

                    if (!this.s[i * n + j]) continue;

                    let u = .5 * (this.u[i * n + j] + this.u[(i + 1) * n + j]);
                    let v = .5 * (this.v[i * n + j] + this.v[i * n + j + 1]);
                    let x = i * h + h2 - dt * u;
                    let y = j * h + h2 - dt * v;
                    this.newM[i * n + j] = this.sampleField(x, y, 'S_FIELD');
                }
            }
            this.m.set(this.newM);
        }

        simulate(dt, grav, iterations) {
            // Reset solid cells
            this.s.fill(0);
            
            // Mark fluid cells
            for(let i = 1; i < this.nx-1; i++) {
                for(let j = 1; j < this.ny-1; j++) {
                this.s[i * this.ny + j] = 1;
                }
            }
            
            for(let obj of objects) {
                if(obj.type === 'svg') {
                    // Existing SVG processing code
                } else { // Handle circle/square
                    const cellsX = Math.floor(obj.x / this.h);
                    const cellsY = Math.floor(obj.y / this.h);
                    const radius = Math.floor(obj.size / (2 * this.h));
                    
                    for(let dx = -radius; dx <= radius; dx++) {
                        for(let dy = -radius; dy <= radius; dy++) {
                            const x = cellsX + dx;
                            const y = cellsY + dy;
                            if(x >= 1 && x < this.nx-1 && y >= 1 && y < this.ny-1) {
                                if(obj.type === 'circle' && dist(dx, dy, 0, 0) > radius) continue;
                                this.s[x * this.ny + y] = 0;
                            }
                        }
                    }
                }
            }
        
            // Add object velocity to fluid
            for(let obj of objects) {
                const velX = obj.x - obj.prevX;
                const velY = obj.y - obj.prevY;
                obj.prevX = obj.x;
                obj.prevY = obj.y;
                this.applyObjectVelocity(obj, velX*10, velY*10);
            }
            // Add pulsating spout
            let pulsation = Math.sin(frameCount * 0.1) * 5;
            let spoutVelocity = this.spoutIntensity + pulsation;
            
            // Spout dimensions (middle 1/3 of left edge)
            let startJ = Math.floor(2*this.ny / 12);
            let endJ = Math.floor( 2*this.ny / 6);
            console.log(startJ, endJ);
            // Inject fluid and velocity
            for(let j = startJ; j < endJ; j++) {
                let index = 1 * this.ny + j; // Leftmost fluid column
                this.m[index] = 1;        // Max density
                this.u[index] = spoutVelocity; // Rightward velocity
            }
    
            this.integrate(dt, grav);
            this.p.fill(0);
            this.uncompress(iterations, dt);
            this.extrapolate();
            this.advectVel(dt);
            this.advectSmoke(dt);
            this.UpdateCells();
        }

        UpdateCells() {
            for(let x = 0; x < this.cellList.length; x++) {
                for(let y = 0; y < this.cellList[x].length; y++) {
                    // Map mass to opacity with clamping
                    let massValue = constrain(this.m[x * this.ny + y], 0, 1);
                    this.cellList[x][y].value = massValue * 255;
                    
                    // Add velocity visualization (optional)
                    let vel = Math.hypot(this.u[x * this.ny + y], this.v[x * this.ny + y]);
                    this.cellList[x][y].color = [0, 50, 200 + vel * 10, massValue * 255];
                }
            }
        }
        applyObjectVelocity(obj, vx, vy) {
            const h = this.h;
            
            let n = this.ny;
            for (let i = 1; i < this.nx-2; i++) {
                for (let j = 1; j < this.ny-2; j++) {
                    // stay only if i, j in object
                    // let dx = i*h + h/2 - obj.x;
                    // let dy = j*h + h/2 - obj.y;
                    // if (Math.abs(dx) >= radius || Math.abs(dy) >= radius) continue;
                    // if (dist(0, 0, dx, dy) >= radius) continue;
                    if (!obj.contains(i*h+h/2, j*h+h/2)) continue;

                    this.s[i*n+j] = 0.;
                    this.m[i*n+j] = 0.;  // try adjusting this val to change trail color
                    this.u[i*n+j] = this.u[(i+1)*n+j] = vx;
                    this.v[i*n+j] = this.v[i*n+j+1] = vy;
                }
            }
        }
    }

class DraggableObject {
    constructor(x, y, type, data) {
        this.x = x;
        this.y = y;
        this.type = type;
        this.data = data;
        this.dragging = false;
        this.prevX = x;
        this.prevY = y;
        
        // Set size based on type
        if(type === 'svg') {
            this.size = max(this.data.w, this.data.h) * 2; // Scale up SVG objects
        } else {
            this.size = 60; // Default size for circles/squares
        }
    }

    show() {
        stroke(255);
        fill(200, 200); // More visible fill
        
        if (this.type === 'svg') {
            push();
            translate(this.x, this.y); // Center position
            scale(0.5); // Match the scale used in parseSVG
            beginShape();
            this.data.points.forEach(p => vertex(p.x, p.y));
            endShape(CLOSE);
            pop();
        } else if (this.type === 'circle') {
            ellipse(this.x, this.y, this.size, this.size);
        } else { // square
            rectMode(CENTER);
            rect(this.x, this.y, this.size, this.size);
        }
    }
    
    contains(px, py) {
        if (this.type === 'svg') {
            const localX = px - (this.x - this.data.w/2);
            const localY = py - (this.y - this.data.h/2);
            return this.data.points.some(p => 
              dist(localX, localY, p.x + this.x - this.data.w/2, p.y + this.y - this.data.h/2) < 5
            );
          }
        if (this.type === 'circle') {
        return dist(px, py, this.x, this.y) < this.size/2;
        } else { // square
        return px > this.x - this.size/2 && px < this.x + this.size/2 &&
                py > this.y - this.size/2 && py < this.y + this.size/2;
        }
    }
    
    update() {
        if (this.dragging) {
        this.x = mouseX;
        this.y = mouseY;
        }
    }
}