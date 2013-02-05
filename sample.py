import random
from Polygon import Polygon
from math import sqrt

def tri_sample(tri):
    trix = [p[0] for p in tri]
    triy = [p[1] for p in tri]
    u = random.uniform(0, 1)
    v = random.uniform(0, 1)
    tmp = sqrt(u);
    a = 1 - tmp;
    b = v * tmp;
    c = 1 - a - b
    sx = a * trix[0] + b * trix[1] + c * trix[2]
    sy = a * triy[0] + b * triy[1] + c * triy[2]
    return (sx, sy)

sample_count = 0
def sample(poly):
    global sample_count
    sample_count += 1 
    a = poly.area() * random.uniform(0, 1)
    strips = poly.triStrip()
    for idx, strip in enumerate(strips):
        for t0 in range(len(strip) - 2):
            tri_area = Polygon(strip[t0:t0+3]).area()
            if tri_area >= a:
                return tri_sample(strip[t0:t0+3])
            a -= tri_area
    # roundoff
    return tri_sample(strips[-1][-3:])

if __name__ == '__main__':
    from pylab import fill, show, plot
    from Polygon.Shapes import Circle
    p = Polygon([(0, 0), (0, 1), (1,1), (1, 0)])
    c = Circle(0.25, (0.3, 0.7))
    sub = p - c
    sub.simplify()
    for strip in sub.triStrip():
        for t0 in range(len(strip) - 2):
            x = [p[0] for p in strip[t0:t0+3]]
            y = [p[1] for p in strip[t0:t0+3]]
            fill(x, y, fc='w')
    samps = [sample(sub) for n in range(10000)]
    x = [p[0] for p in samps]
    y = [p[1] for p in samps]
    plot(x, y, 'k.')
    show()
