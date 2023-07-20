                                                            /*
JSRENDER é um renderizador puramente implementado em software
usando poucas diretivas canvas do HTML5.
G.G.
*/
class Mapping
{
    constructor(Xwmin, Ywmin, Xwmax, Ywmax, Xsmin, Ysmin, Xsmax, Ysmax, isotropic=false)
    {
        this.Xwmin = Xwmin;
        this.Xwmax = Xwmax;
        this.Ywmin = Ywmin;
        this.Ywmax = Ywmax;

        this.Xsmin = Xsmin;
        this.Xsmax = Xsmax;
        this.Ysmin = Ysmin;
        this.Ysmax = Ysmax;
        /*
            larguraPixel = larguraCena/numeroPixels;

            numeroPixelsNaHorizontal = larguraCena/larguraPixel;
            numeroPixelsNaHorizontal = larguraRegiao/larguraPixel;

            O mesmo vale para a aultura.

            O pixel pode ser anisotrópico, quando larguraPixel != alturaPixel
            ou isotrópico quando larguraPixel = alturaPixel.

            Se quisermos forçar a isotropia,
            devemos pegar o maior tamanho de pixel, pois
            assim estaremos garantido pegar a menor quantidade de pixels
            considerando a altura e a largura.
        */
        this.pixelWidth = (Xwmax-Xwmin)/(Xsmax-Xsmin);
        this.pixelHeight = (Ywmax - Ywmin)/(Ysmax - Ysmin);
        
        if  (isotropic)
        {
            var ms = Math.max(this.pixelWidth, this.pixelHeight);
            this.pixelWidth = ms;
            this.pixelHeight = ms;
        }
    }

    xToSc(x)
    {
        return this.Xsmin + (x-this.Xwmin)/this.pixelWidth;
    }
    
    yToSc(y)
    {
        return this.Ysmin + ( this.Ysmax -  ((y - this.Ywmin)/this.pixelHeight) );
    }

    ixToWindow(ix)
    { //ix = this.Xsmin + (x-this.Xwmin)/this.pixelWidth;
      // x = (ix - Xsmin) * pW + Xwin  
        return (ix - this.Xsmin) * this.pixelWidth + this.Xwmin; 
    }

    iyToWindow(iy)
    {//iy = this.Ysmin + ( this.Ysmax -  ((y - this.Ywmin)/this.pixelHeight) )
        //y = ( (iy + this.Ysmin) + this.Ysmax ) * pW - this.Ywmin
        return ((this.Ysmin - iy) + this.Ysmax) * this.pixelHeight + this.Ywmin; 
    }

    drawBorder(context, color="red")
    {
        var st = context.style;
        context.beginPath()
        context.strokeStyle = color;
        context.moveTo(this.Xsmin, this.Ysmin);
        context.lineTo(this.Xsmax, this.Ysmin);
        context.lineTo(this.Xsmax, this.Ysmax);
        context.lineTo(this.Xsmin, this.Ysmax);
        context.lineTo(this.Xsmin, this.Ysmin);
        context.stroke();
        context.strokeStyle = st;
    }

    drawLine(context, x1, y1, x2, y2, color="black")
    {
        //console.log(this.xToSc(x1) + ", " + this.yToSc(y1) + " ::: " + this.xToSc(x2) + ",  " + this.yToSc(y2));
        var st = context.style;
        context.beginPath();
        context.strokeStyle = color;
        context.moveTo(this.xToSc(x1), this.yToSc(y1));
        context.lineTo(this.xToSc(x2), this.yToSc(y2));
        context.stroke();
        context.strokeStyle = st;
    }
}

function desenharCirculo(context, m, c, raio, amostras=100)
{
    context.beginPath();
    var delta = (2 * Math.PI)/amostras;
    var ang = 0;
    var fx = 0;
    var fy = 0;
    for (var i = 0; i < amostras; i++)
    {
        var xi = c[0] + raio * Math.cos(ang);
        var yi = c[1] + raio * Math.sin(ang);
        ang += delta;
        if (i == 0)
        {
            context.moveTo(m.xToSc(xi), m.yToSc(yi));
            fx = xi;
            fy = yi;
        }
        else
        {
            context.lineTo(m.xToSc(xi), m.yToSc(yi));
        }
    }
    context.lineTo(m.xToSc(fx), m.yToSc(fy));
    context.stroke();
}

class Vector2D 
{
    constructor(x = 0, y = 0)
    {
        this.x = x;
        this.y = y;
    }

    sum(vector) 
    {
        return new Vector2D(this.x + vector.x, this.y + vector.y);
    }

    sub(vector)
    {
        return new Vector2D(this.x - vector.x, this.x - vector.y);
    }

    dot(v)
    {
        return this.x * v.x + this.y * v.y;
    }

    det(v)
    {
        this.x * v.y - this.y * v.x;
    }
}

class Vector3D 
{
    constructor(x = 0, y = 0, z = 0, w = 0)
    {
        this.x = x;
        this.y = y;
        this.z = z;
        this.w = w;
    }

    sum(vector) 
    {
        return new Vector3D(this.x + vector.x, this.y + vector.y, this.z + vector.z);
    }

    sub(vector)
    {
        return new Vector3D(this.x - vector.x, this.y - vector.y, this.z - vector.z);
    }

    dot(v)
    {
        return this.x * v.x + this.y * v.y + this.z * v.z;
    }

    scale(k)
    {
        return new Vector3D(k * this.x, k * this.y, k * this.z, 0);
    }

    size()
    {
        return Math.sqrt(this.x * this.x + this.y * this.y + this.z * this.z);
    }

    cross(v)
    {
        var x = this.y * v.z - this.z * v.y;
        var y = this.z * v.x - this.x * v.z;
        var z = this.x * v.y - this.y * v.x;
        return new Vector3D(x, y, z); 
    }

    normalize()
    {
        var n = this.size();
        if (n != 0)
        {
            return new Vector3D(this.x/n, this.y/n, this.z/n, this.w);
        }
        else
        {
            return new Vector3D(this.x, this.y, this.z, this.w);
        }
    }
}

class Point2D
{
    constructor(x = 0, y = 0)
    {
        this.x = x;
        this.y = y;       
    }

    sum(vector) 
    {
        return new Point3D(this.x + vector.x, this.y + vector.y);
    }

    sub(point)
    {
        return new Vector2D(this.x - point.x, this.y - point.y);
    }
}

class Point3D
{
    constructor(x = 0, y = 0, z = 0, w = 1)
    {
        this.x = x;
        this.y = y;
        this.z = z;
        this.w = w;
    }

    sum(vector) 
    {
        return new Point3D(this.x/this.w + vector.x, this.y/this.w + vector.y, this.z/this.w + vector.z, 1);
    }

    sub(point)
    {
        var w = point.w;
        return new Vector3D(this.x/this.w - point.x/w, this.y/this.w - point.y/w, 
                            this.z/this.w - point.z/w, 0);
    }
}

class Line2
{
    constructor(p1, p2, color='black')
    {
        this.begin = p1;
        this.end = p2;
        this.color = color;
    }
}

class Matrix2
{
    constructor(data = null)
    {
        if (data == null)
        {
            this.data = [[1, 0], [0, 1]];
        }
        else
        {
            this.data = [[0, 0], [0, 0]];
            for (var i = 0; i < data.length; i++)
            {
                for (var j = 0; j < data[i].length; j++)
                {
                    this.data[i][j] = data[i][j]; 
                }
            }    
        }
    }

    multP(p)
    {
        return new Point2D(this.data[0][0] * p.x + this.data[0][1] * p.y,
                            this.data[1][0] * p.x + this.data[1][1] * p.y);  
    }

    multV(v)
    {
        return new Vector2D(this.data[0][0] * v.x + this.data[0][1] * v.y,
            this.data[1][0] * v.x + this.data[1][1] * v.y);  
    }

    mult(m)
    {
        return new Matrix2([ [this.data[0][0] * m.data[0][0] + this.data[0][1] * m.data[1][0],
                           this.data[0][0] * m.data[0][1] + this.data[0][1] * m.data[1][1]],
                           [this.data[1][0] * m.data[0][0] + this.data[1][1] * m.data[1][0],
                           this.data[1][0] * m.data[0][1] + this.data[1][1] * m.data[1][1]] ])
    }
}

class Matrix3
{
    constructor(data = null)
    {
        if (data == null)
        {
            this.data = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
        }
        else
        {
            this.data = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];
            for (var i = 0; i < data.length; i++)
            {
                for (var j = 0; j < data[i].length; j++)
                {
                    this.data[i][j] = data[i][j]; 
                }
            }    
        }
    }

    multP(p)
    {
        return new Point3D(this.data[0][0] * p.x + this.data[0][1] * p.y + this.data[0][2] * p.z,
                            this.data[1][0] * p.x + this.data[1][1] * p.y + this.data[1][2] * p.z,
                            this.data[2][0] * p.x + this.data[2][1] * p.y + this.data[2][2] * p.z,);  
    }

    multV(v)
    {
        return new Vector3D(this.data[0][0] * v.x + this.data[0][1] * v.y + this.data[0][2] * v.z,
            this.data[1][0] * v.x + this.data[1][1] * v.y + this.data[1][2] * v.z,
            this.data[2][0] * v.x + this.data[2][1] * v.y + this.data[2][2] * v.z);  
    }

    mult(m)
    {
        let a = this.data;
        let b = m.data;
        return new Matrix3([[a[0][0] * b[0][0] + a[0][1] * b[1][0] + a[0][2] * b[2][0],
                           a[0][0] * b[0][1] + a[0][1] * b[1][1] + a[0][2] * b[2][1],
                           a[0][0] * b[0][2] + a[0][1] * b[1][2] + a[0][2] * b[2][2]],
                           [a[1][0] * b[0][0] + a[1][1] * b[1][0] + a[1][2] * b[2][0],
                           a[1][0] * b[0][1] + a[1][1] * b[1][1] + a[1][2] * b[2][1],
                           a[1][0] * b[0][2] + a[1][1] * b[1][2] + a[1][2] * b[2][2]],
                           [a[2][0] * b[0][0] + a[2][1] * b[1][0] + a[2][2] * b[2][0],
                           a[2][0] * b[0][1] + a[2][1] * b[1][1] + a[2][2] * b[2][1],
                           a[2][0] * b[0][2] + a[2][1] * b[1][2] + a[2][2] * b[2][2]]]
                           )
    }
}

class Circle
{
    constructor(x, y, r)
    {
        this.x = x;
        this.y = y;
        this.radio = r;
        this.origin = new Point3D(0, 0, 1);
        this.model = new Matrix3([[1, 0, 0], [0, 1, 0], [0, 0, 1]]);
        this.resolution = 3;
    }

    render()
    {
        let delta = (2*Math.PI)/(this.resolution);
        let ang = delta;
        let x = this.origin.x +  Math.cos(ang) * this.radio;
        let y = this.origin.y +  Math.sin(ang) * this.radio;
        let first = new Point2D(x, y);
        let begin = new Point2D(x, y);
        let lines = []
        for (var i = 1; i < this.resolution; i++)
        {
            let x = this.origin.x +  Math.cos(ang) * this.radio;
            let y = this.origin.y + Math.sin(ang) * this.radio;
            ang += delta;
            let end = new Point3D(x, y, 1);
            lines.push(new Line2(begin, end));
            begin = end;
        }
        let end = first;
        lines.push(new Line2(end, begin));
        return lines;
    }
}

class Triangle
{
    constructor(x1, y1, x2, y2, x3, y3)
    {
        this.p1 = new Point3D(x1, y1, 1);
        this.p2 = new Point3D(x2, y2, 1);
        this.p3 = new Point3D(x3, y3, 1);
        this.origin = new Point3D(0, 0, 1);
        this.model = new Matrix3([[1, 0, 0], [0, 1, 0], [0, 0, 1]]);
    }

    render()
    {
        return [new Line2(this.p1, this.p2), new Line2(this.p2, this.p3), new Line2(this.p3, this.p1)];
    }
}

class Polygon
{
    constructor(xs, ys)
    {
        this.points = []
        this.n = Math.min(xs.length, ys.length);
        for (var i = 0; i < this.n; i++)
        {
            this.points.push(new Point3D(xs[i], ys[i], 1))
        }
        this.origin = new Point3D(0, 0, 1);
        this.model = new Matrix3([[1, 0, 0], [0, 1, 0], [0, 0, 1]]);
        
        this.lines = []
        for (var j = 0; j < this.n-1; j++)
        {
            this.lines.push(new Line2(this.points[j], this.points[j+1]));
        }
        this.lines.push(new Line2(this.points[this.n-1], this.points[0]));
    }

    render()
    {
        return this.lines;
    }
}

class Window
{
    constructor(ctx, xmin, ymin, xmax, ymax)
    {
        this.xmin = xmin;
        this.ymin = ymin;
        this.xmax = xmax;
        this.ymax = ymax;
        this.origin = new Point3D(0, 0, 0);
        this.mapping = new Mapping(xmin, ymin, xmax, ymax, 0, 0, ctx.canvas.width, ctx.canvas.height, true);
        //this.orientation =  new Matrix2([[0.707, -0.707], [0.707, 0.707]]);
        this.orientation =  new Matrix3([[1, 0, 0], [0, 1, 0], [0, 0, 1]]); 
        this.objects = [];
        this.context = ctx;
    }

    toWindow(p, obj)
    {
        var tp1 = obj.origin
        var tp2 = this.orientation.mult(obj.model).multP(p);
        return tp1.sum(tp2);
    }

    draw()
    {
        for (var i = 0; i < this.objects.length; i++)
        {
            let obj = this.objects[i];
            let lines = obj.render();
            for (var j = 0; j < lines.length; j++)
            {
                let line = lines[j]
                let p1 = this.toWindow(line.begin, obj);
                let p2 = this.toWindow(line.end, obj);
                this.mapping.drawLine(this.context, p1.x, p1.y, p2.x, p2.y);
            }
        }
        this.mapping.drawBorder(this.context);
    }
}

class Window3D
{
    constructor(ctx, xmin, ymin, xmax, ymax, camera, projection)
    {
        this.xmin = xmin;
        this.ymin = ymin;
        this.xmax = xmax;
        this.ymax = ymax;
        this.origin = new Point3D(0, 0, 0);
        this.mapping = new Mapping(xmin, ymin, xmax, ymax, 0, 0, ctx.canvas.width, ctx.canvas.height, true);
        this.orientation =  new Matrix4([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]]); 
        this.objects = [];
        this.context = ctx;
        this.camera = camera
        this.projection = projection
    }

    toWindow(p, obj)
    {
        var tp2 = this.projection.mult(this.camera.matrix).mult(obj.model).mult(this.orientation).multP(p);
        let x = tp2.x/tp2.w;
        let y = tp2.y/tp2.w;
        let z = tp2.z/tp2.w;
        console.log(x + ", " + y + ", " + z);
        let np = new Point3D(x, y, z);
        return np;
    }

    draw()
    {
        for (var i = 0; i < this.objects.length; i++)
        {
            let obj = this.objects[i];
            let lines = obj.render();
            for (var j = 0; j < lines.length; j++)
            {
                let line = lines[j]
                let p1 = this.toWindow(line.begin, obj);
                let p2 = this.toWindow(line.end, obj);
                if (outofrange(p1) || outofrange(p2))
                {
                    continue;
                }
                this.mapping.drawLine(this.context, p1.x, p1.y, p2.x, p2.y, line.color);
            }
        }
        this.mapping.drawBorder(this.context);
    }
}

function outofrange(p)
{
    return (p.x > 1 || p.x < -1) || (p.y > 1 || p.y < -1) || (p.z > 1 || p.z < -1); 
}

function rotate(theta)
{
    r = theta * Math.PI/180.0
    c = Math.cos(r)
    s = Math.sin(r)
    return new Matrix3([[c, -s, 0], [s, c, 0], [0, 0, 1]])
}

function rotate3x(theta)
{
    r = theta * Math.PI/180.0
    c = Math.cos(r)
    s = Math.sin(r)
    return new Matrix4([[1, 0, 0, 0], [0, c, -s, 0], [0, s, c, 0], [0, 0, 0, 1]]);
}

function rotate3y(theta)
{
    r = theta * Math.PI/180.0
    c = Math.cos(r)
    s = Math.sin(r)
    return new Matrix4([[c, 0, s, 0], [0, 1, 0, 0], [-s, 0, c, 0], [0, 0, 0, 1]]);
}

function rotate3z(theta)
{
    r = theta * Math.PI/180.0
    c = Math.cos(r)
    s = Math.sin(r)
    return new Matrix4([[c, -s, 0, 0], [s, c, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]]);
}

function translate(dx, dy)
{
    return new Matrix3([[1, 0, dx], [0, 1, dy], [0, 0, 1]])
}

function translate3(dx, dy, dz)
{
    return new Matrix4([[1, 0, 0, dx], [0, 1, 0, dy], [0, 0, 1, dz], [0, 0, 0, 1]])
}

function scale(sx, sy)
{
    return new Matrix3([[sx, 0, 0], [0, sy, 0], [0, 0, 1]])
}

function shear(sx, sy)
{
    return new Matrix3([[1, sx, 0], [sy, 1, 0], [0, 0, 1]])
}
4
function area2(a, b, c)
{
    return (a.x - c.x) * (b.y - c.y) - 
           (a.y - c.y) * (b.x - c.x);
}

function  ccw(points)
{
    let k = 0;
    for (var i = 1; i < points.length; i++)
    {
        if ( 
            (points[i].x <= points[k].x) &&
            (points[i].x < points[k].x || points[i].y  < p[k].y) 
            )
        {
            k = i;
        }
    }
    let prev = k - 1, next = k + 1;
    if (prev < 0) prev = points.length - 1;
    if (next >= points.length) next = 0;
    return area2(points[prev], points[k], points[next])
}

function insideTriangle(t, p, mywindow, eps=1.0e-18)
{
    p1 = mywindow.toWindow(t.p1, t);
    p2 = mywindow.toWindow(t.p2, t);
    p3 = mywindow.toWindow(t.p3, t);
    return area2(p1, p2, p) >= eps &&
           area2(p2, p3, p) >= eps &&
           area2(p3, p1, p) >= eps;
}

function insidePolygon(polygon, p, mywindow, eps=1.0e-18)
{
    var edges = polygon.render()
    let counter = 0;
    for (var i = 0; i < edges.length; i++)
    {
        e = edges[i];

        let a = mywindow.toWindow(e.begin, polygon);
        let b = mywindow.toWindow(e.end, polygon); 
        if ( 
            (p.y <= b.y && p.y > a.y && area2(a, b, p) > 0) || 
                (p.y <= a.y && p.y > b.y && area2(b, a, p) > 0) 
            )
        {
            counter++;
        }
    }
    return counter % 2 != 0 && counter > 0;
}

class Matrix4
{
    constructor(data = null)
    {
        if (data == null)
        {
            this.data = [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]];
        }
        else
        {
            this.data = [[0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]];
            for (var i = 0; i < data.length; i++)
            {
                for (var j = 0; j < data[i].length; j++)
                {
                    this.data[i][j] = data[i][j]; 
                }
            }
        }
    }

    multP(p)
    {
        return new Point3D(this.data[0][0] * p.x + this.data[0][1] * p.y + this.data[0][2] * p.z + this.data[0][3] * p.w,
                            this.data[1][0] * p.x + this.data[1][1] * p.y + this.data[1][2] * p.z + this.data[1][3] * p.w,
                            this.data[2][0] * p.x + this.data[2][1] * p.y + this.data[2][2] * p.z + this.data[2][3] * p.w,
                            this.data[3][0] * p.x + this.data[3][1] * p.y + this.data[3][2] * p.z + this.data[3][3] * p.w);  
    }

    multV(p)
    {
        return new Vector3D(this.data[0][0] * p.x + this.data[0][1] * p.y + this.data[0][2] * p.z + this.data[0][3] * p.w,
            this.data[1][0] * p.x + this.data[1][1] * p.y + this.data[1][2] * p.z + this.data[1][3] * p.w,
            this.data[2][0] * p.x + this.data[2][1] * p.y + this.data[2][2] * p.z + this.data[2][3] * p.w,
            this.data[3][0] * p.x + this.data[3][1] * p.y + this.data[3][2] * p.z + this.data[3][3] * p.w);   
    }

    mult(m)
    {
        let a = this.data;
        let b = m.data;
        return new Matrix4([[a[0][0] * b[0][0] + a[0][1] * b[1][0] + a[0][2] * b[2][0] + a[0][3] * b[3][0],
                           a[0][0] * b[0][1] + a[0][1] * b[1][1] + a[0][2] * b[2][1] + a[0][3] * b[3][1],
                           a[0][0] * b[0][2] + a[0][1] * b[1][2] + a[0][2] * b[2][2] + a[0][3] * b[3][2],
                           a[0][0] * b[0][3] + a[0][1] * b[1][3] + a[0][2] * b[2][3] + a[0][3] * b[3][3]],
                           [a[1][0] * b[0][0] + a[1][1] * b[1][0] + a[1][2] * b[2][0] + a[1][3] * b[3][0],
                           a[1][0] * b[0][1] + a[1][1] * b[1][1] + a[1][2] * b[2][1] + a[1][3] * b[3][1],
                           a[1][0] * b[0][2] + a[1][1] * b[1][2] + a[1][2] * b[2][2] + a[1][3] * b[3][2],
                           a[1][0] * b[0][3] + a[1][1] * b[1][3] + a[1][2] * b[2][3] + a[1][3] * b[3][3]],
                           [a[2][0] * b[0][0] + a[2][1] * b[1][0] + a[2][2] * b[2][0] + a[2][3] * b[3][0],
                           a[2][0] * b[0][1] + a[2][1] * b[1][1] + a[2][2] * b[2][1] + a[2][3] * b[3][1],
                           a[2][0] * b[0][2] + a[2][1] * b[1][2] + a[2][2] * b[2][2] + a[2][3] * b[3][2],
                           a[2][0] * b[0][3] + a[2][1] * b[1][3] + a[2][2] * b[2][3] + a[2][3] * b[3][3]],
                           [a[3][0] * b[0][0] + a[3][1] * b[1][0] + a[3][2] * b[2][0] + a[3][3] * b[3][0],
                           a[3][0] * b[0][1] + a[3][1] * b[1][1] + a[3][2] * b[2][1] + a[3][3] * b[3][1],
                           a[3][0] * b[0][2] + a[3][1] * b[1][2] + a[3][2] * b[2][2] + a[3][3] * b[3][2],
                           a[3][0] * b[0][3] + a[3][1] * b[1][3] + a[3][2] * b[2][3] + a[3][3] * b[3][3]]]
                           )
    }
}

class Camera
{
    constructor(position, viewpoint, up)
    {
        this.position = position;
        this.n = position.sub(viewpoint).normalize();
        this.v = up.sub( this.n.scale( this.n.dot(up) ) ).normalize();
        this.u = this.v.cross(this.n).normalize();

        var R = new Matrix4( [ [this.u.x, this.u.y, this.u.z, 0], [this.v.x, this.v.y, this.v.z, 0], [this.n.x, this.n.y, this.n.z, 0], [0, 0, 0, 1]  ] );
        var T = new Matrix4( [ [1, 0, 0, -this.position.x], [0, 1, 0, -this.position.y], [0, 0, 1, -this.position.z], [0, 0, 0, 1]] )
        
        this.matrix = R.mult(T);
    }
}

class Cube
{
    constructor(center, size=1)
    {
        this.center = center;
        this.size = size;
        this.lines = [];
        this.model = new Matrix4();
        this.update();
    }

    update()
    {
        let xmin = this.center.x - this.size/2;
        let xmax = this.center.x + this.size/2;
        let ymin = this.center.y - this.size/2;
        let ymax = this.center.y + this.size/2;
        let zmin = this.center.z - this.size/2;
        let zmax = this.center.z + this.size/2;

        let p1 = new Point3D(xmin, ymin, zmin); //BACK, LEFT, BOT
        let p2 = new Point3D(xmax, ymin, zmin); //BACK, RIGHT, BOT
        let p3 = new Point3D(xmax, ymax, zmin); //BACK, RIGHT, TOP
        let p4 = new Point3D(xmin, ymax, zmin); //BACK, LEFT, TOP
        
        let p5 = new Point3D(xmax, ymin, zmax); //FRONT, RIGHT, BOT
        let p6 = new Point3D(xmax, ymax, zmax); //FRONT, RIGHT, TOP
        let p7 = new Point3D(xmin, ymax, zmax); //FRONT, LEFT, TOP
        let p8 = new Point3D(xmin, ymin, zmax); //FRONT, LEFT, BOT
        
        //BOTTOM
        var cb = "yellow";
        this.lines.push(new Line2(p5, p8, cb));
        this.lines.push(new Line2(p8, p1, cb));
        this.lines.push(new Line2(p1, p2, cb));
        this.lines.push(new Line2(p2, p5, cb));

        //RIGHT
        var cr = "blue";
        this.lines.push(new Line2(p5, p6, cr));
        this.lines.push(new Line2(p6, p3, cr));
        this.lines.push(new Line2(p3, p2, cr));
        this.lines.push(new Line2(p2, p5, cr));

        //LEFT
        var cl = "red";
        this.lines.push(new Line2(p8, p7, cl));
        this.lines.push(new Line2(p7, p4, cl));
        this.lines.push(new Line2(p4, p1, cl));
        this.lines.push(new Line2(p1, p8, cl));

        //FRONT FACE
        var cf = "black";
        this.lines.push(new Line2(p5, p6, cf));
        this.lines.push(new Line2(p6, p7, cf));
        this.lines.push(new Line2(p7, p8, cf));
        this.lines.push(new Line2(p8, p5, cf));

        //BACK FACE
        var ca = "cyan";
        this.lines.push(new Line2(p1, p2, ca));
        this.lines.push(new Line2(p2, p3, ca));
        this.lines.push(new Line2(p3, p4, ca));
        this.lines.push(new Line2(p4, p1, ca));

        //TOP FACE
        var ct = "magenta";
        this.lines.push(new Line2(p3, p4, ct));
        this.lines.push(new Line2(p4, p7, ct));
        this.lines.push(new Line2(p7, p6, ct));
        this.lines.push(new Line2(p6, p3, ct));
    }

    render()
    {
        return this.lines;
    }
}