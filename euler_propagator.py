def magnitude(v):
    return (v[0]**2 + v[1]**2 + v[2]**2)**0.5

def unit_vector(v):
    return [i/magnitude(v) for i in v]

pos = [3207,5459,2714]
vel = [-6.532,0.7835,6.142]
r = magnitude(pos)

r_prev = 0

t = 0
T = 1 # time step

while r_prev < r:
    r_prev = r
    pos[0] += vel[0]*T
    pos[1] += vel[1]*T
    pos[2] += vel[2]*T

    a = 3.986e5 / r**2
    acc = [-i*a for i in unit_vector(pos)]

    vel[0] += acc[0]*T
    vel[1] += acc[1]*T
    vel[2] += acc[2]*T

    r = magnitude(pos)
    t += 1

    print(f'{t},{r},{vel},{pos}')
