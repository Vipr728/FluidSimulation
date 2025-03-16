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
        const svgFiles = ['abhi.txt', 'cow.txt', 'space.txt', 'george.txt', 'arrow.txt'];
        svgFiles.forEach(filename => {
            // Use loadXML instead of loadStrings
            svgTemplates[filename] = loadStrings(`svgs/${filename}`);
            console.log(svgTemplates[filename]); // Log loaded SVGs
        
        });
    }
    
    function addSVGObject(filename) {
        const svg = svgTemplates[filename];
        let svgString = svg;
        console.log("svgString" + svgString); // Log the SVG string
        const svgData = parseSVG(svgString, 0, 0);
        console.log("SVG data" + svgData); // Log the parsed SVG data
        objects.push(new DraggableObject(mouseX, mouseY, 'svg', svgData));
    }

    function setup() {
        createCanvas(800, 520);
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
        let spoutSlider = createSlider(0,1000, 250, 1);
        spoutSlider.input(() => {
            gridInstance.spoutIntensity = spoutSlider.value();
        });
        let addCircle = createButton('Add Circle');
        addCircle.mousePressed(() => addNewObject('circle'));
        let addSquare = createButton('Add Square');
        addSquare.mousePressed(() => addNewObject('square'));
        let addArc = createButton('Add SemiCircle');
        addArc.mousePressed(() => addNewObject('arc'));

        function addNewObject(type) {
            objects.push(new DraggableObject(width/2, height/2, type, 60));
        }
        Object.keys(svgTemplates).forEach(name => {
            let btn = createButton(`Add ${name}`);
            btn.mousePressed(() => addSVGObject(name));}
        );

        objects.push(new DraggableObject(width/2, height/2, 'png', {
            w: 80,
            h: 80,
            path: '/FluidSim/pngs/flash.png',
            mass: 2,
        }))
        // objects.push(new DraggableObject(width/2, height/2, 'png', {
        //     w: 80,
        //     h: 80,
        //     path: '/FluidSim/pngs/cup.png',
        //     mass: 3,
        // }))
    }

    function draw() {
        background('#633f33');
        // Draw fluid first
        drawCells(gridInstance.cellList);
        gridInstance.simulate(0.1, 9.81, 20);
        
        // Then draw objects on top
        for(let obj of objects) {
            obj.update();
            obj.show();
        }

        // image(img, 80, 40, 20, 20, 0, 0, img.width, img.height)
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
                let velX, velY;
                if (obj.dragging) {
                    velX = (obj.x - obj.prevX);
                    velY = (obj.y - obj.prevY);
                } else {
                    [velX, velY] = this.applyPressureToObject(obj);
                    obj.x += velX;
                    obj.y += velY;
                }
                obj.prevX = obj.x;
                obj.prevY = obj.y;

                this.applyObjectVelocity(obj, velX*10, velY*10);
                obj.calcTorque(this);
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

            let secondSpoutI = this.nx - 1; // Rightmost column (excluding boundary)
            let secondSpoutStartJ = Math.floor(3 * this.ny / 4); // Start at 3/4 of the grid height
            let secondSpoutEndJ = this.ny - 2; // End at the bottom row (excluding boundary)

            // Inject fluid and velocity for the second spout
            for(let j = secondSpoutStartJ; j < secondSpoutEndJ; j++) {
                let index = secondSpoutI * this.ny + j;
                this.m[index] = 1; // Max density

                // Set inward velocity (negative u and v for whirlpool effect)
                this.u[index] = -spoutVelocity * 0.5; // Horizontal velocity (inward)
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
                    // this.cellList[x][y].color = [0, 200, 100+vel*20, massValue * 255];
                    this.cellList[x][y].color = [233, 220, 211, massValue * 255];
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
                    obj.vx = vx;
                    obj.vy=vy;
                }
            }
        }
        applyPressureToObject(obj) {
            const h = this.h;
            let n = this.ny;

            let velX = 0, velY = 0;

            for (let i = 1; i < this.nx-2; i++) {
                for (let j = 1; j < this.ny-2; j++) {
                    if (!obj.contains(i*h+h/2, j*h+h/2)) continue;

                    for (let [dx, dy] of [[-1, 0], [1, 0], [0, 1], [0, -1]]) {
                        let i2 = i+dx, j2 = j+dy;
                        if (this.s[i2*n+j2] == 0.) continue;
                        // move by -dx, -dy
                        let mag = this.p[i2*n+j2]/(320000000*obj.mass);
                        velX -= dx*mag;
                        velY -= dy*mag;
                    }
                }
            }
            return [velX, velY]
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
        this.drag = 0;
        this.vx=0;
        this.vy=0;
        this.angle=0;
        this.angularVelocity=0;
        this.torque = 0;
        this.mass = 1.0;

        if (this.data.mass != null) this.mass = this.data.mass;
        
        if(type === 'svg') {
            this.size = max(this.data.w, this.data.h) * 2; // Scale up SVG objects
        } else if (type == 'png') {
            this.img = loadImage(this.data.path);
        } else {
            this.size = 60;
        }
    }

    show() {
        stroke(255);
        fill(225, 225); // More visible fill
        push();
        translate(this.x, this.y); // Move to object's position
        rotate(this.angle); // Apply rotation
        if (this.type === 'svg') {
            scale(0.5); // Match the scale used in parseSVG
            beginShape();
            this.data.points.forEach(p => vertex(p.x, p.y));
            endShape(CLOSE);
        } else if (this.type == 'png') {
            console.log(`x=${this.x} y=${this.y} w=${this.data.w} h=${this.data.h} width=${width} height=${height}`)
            image(this.img, -this.data.w/2, -this.data.h/2, this.data.w, this.data.h, 0, 0, this.img.width, this.img.height)
        } else if (this.type === 'circle') {
            ellipse(0, 0, this.size, this.size); // Draw at (0, 0) after translate
        } else if (this.type === 'arc') {
            arc(0, 0, this.size, this.size, 0, PI);
        } else { // square
            rectMode(CENTER);
            rect(0, 0, this.size, this.size); // Draw at (0, 0) after translate
        }
        pop();
    }
    
    contains(px, py) {
        px -= this.x;
        py -= this.y;

        let a = 0;  // this.angle
        let c = Math.cos(a), s = Math.sin(a);
        [px, py] = [px*c-py*s, py*c+px*s];

        px += this.x;
        py += this.y;

        if (this.type == 'png') {
            let x = this.x-this.data.w/2;
            let y = this.y-this.data.h/2;
            let localX = (px-x)/this.data.w*this.img.width;
            let localY = (py-y)/this.data.h*this.img.height;

            if (localX < 0 || localX >= this.img.width) return false;
            if (localY < 0 || localY >= this.img.height) return false;
            let pixel = this.img.get(localX, localY);
            let alpha = pixel[3];
            return alpha > 2;
        }

        if (this.type === 'svg') {
            const localX = px - (this.x - this.data.w/2);
            const localY = py - (this.y - this.data.h/2);
            return this.data.some(p => 
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

    getDrag() {
        this.drag = 2*0/1;

    }
    
    update() {
        if (this.dragging) {
            this.x = mouseX;
            this.y = mouseY;
        } else {
            // Update position based on velocity
            this.x += this.vx;
            this.y += this.vy;
    
            // Update rotation based on angular velocity
            this.angle += this.angularVelocity;
    
            // Apply damping to angular velocity (to slow down rotation over time)
            this.angularVelocity *= 0.98;
        }
    }
    calcTorque(fluid) {
        const h = fluid.h;
        let n = fluid.ny;
        this.torque = 0;
    
        for (let i = 1; i < fluid.nx - 2; i++) {
            for (let j = 1; j < fluid.ny - 2; j++) {
                if (!this.contains(i * h + h / 2, j * h + h / 2)) continue;
    
                // Calculate the vector from the object's center to the fluid cell
                let dx = (i * h + h / 2) - this.x;
                let dy = (j * h + h / 2) - this.y;
    
                // Calculate the cross product of the position vector and velocity vector
                // Torque = r x F = dx * forceY - dy * forceX
                let forceX = fluid.u[i * n + j];
                let forceY = fluid.v[i * n + j];
                this.torque += dx * forceY - dy * forceX;
            }
        }
    
        // Convert torque to angular acceleration (torque / moment of inertia)
        // For simplicity, assume moment of inertia is proportional to size^2
        let momentOfInertia = this.size * this.size * 0.5 * 1000;
        this.angularVelocity += this.torque / momentOfInertia;
    }
}