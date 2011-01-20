import numpy as np


def r(phi):
        return 0
        
def D(phi, sharpness=1):
    if 0:
        def func(x):
        
            if np.abs(x) >= (np.pi-.001):
                x = np.sign(x)*(np.pi-.001)    
        
            return np.sin(x)

    if 1:
        def func(x):
            if np.abs(x) >= (np.pi-.001):
                x = np.sign(x)*(np.pi-.001)   
            alpha = -1*np.log(1/np.pi)
            val = np.sin( 1 / np.exp( sharpness*np.abs(x) - alpha ) )*np.sign(x)
            return val

    if type(phi) is list:
        val = [func(x) for x in phi]
        return val
    else:
        val = func(phi)
        return val 
    
def N(phi):
    return 0
    
    
class Post:
    def __init__(self, pos=[0,0], radius=1):
        self.pos = np.array(pos)
        self.radius = radius


class Fly:
    def __init__(self, post, initial_pos=[0,0], initial_vel=[0,0], initial_ori=None, dt=0.01, simlen=100):
        self.post = post        
        self.dt = dt
        self.step = 0
        self.simlen = simlen
        self.time = np.arange(0,simlen*dt,dt)
        
        if initial_ori is None: 
            initial_ori = np.arctan2(initial_vel[1], initial_vel[0])
        
        # DYNAMICS
        self.rotmom = 1.5e-2 #kg/m^2 (housefly)
        self.k = 0.18e-1 #kg*m^2/s (housefly) 
        
        # STATES
        self.pos = np.zeros([simlen, 2])
        self.vel = np.zeros([simlen, 2])
        self.world_ori = np.zeros([simlen, 1])
        self.speed = np.zeros([simlen, 1])
        
        # RELATIVE states
        self.angle_to_post = np.zeros([simlen, 1])
        self.dist_to_post = np.zeros([simlen, 1])
        self.angle_subtended = np.zeros([simlen, 1])
        self.expansion = np.zeros([simlen, 1])
        
        # INITIAL CONDITIONS
        self.pos[0] = initial_pos
        self.vel[0] = initial_vel
        self.world_ori[0] = initial_ori
        
        # MODES
        self.legext = False
        self.fixate = True
        self.hover = False
        self.cruise = False
        
        self.calc_relative_states(0)
    
    def calc_relative_states(self, step, dt=None):
        if dt is None:
            dt = self.dt
            
        # ANGLE to post
        vec = self.post.pos-self.pos[step]
        world_angle_to_post = np.arctan2( vec[1], vec[0] )
        angle = world_angle_to_post - self.world_ori[step]
        self.angle_to_post[step] = angle
        self.angle_to_post_diff = np.diff(self.angle_to_post)
            
        # SPEED
        self.speed[step] = np.linalg.norm(self.vel[step])        
                    
        # DIST to post
        self.dist_to_post[step] = np.linalg.norm(vec)
            
        # ANGLE subtended by post
        angle = 2*np.arctan2(self.post.radius, self.dist_to_post[step])
        angle = np.nan_to_num(angle)
        self.angle_subtended[step] = angle
            
        # EXPANSION
        if self.angle_subtended[step-1] != 0:
            optic_speed = (self.angle_subtended[step]-self.angle_subtended[step-1])/dt
            self.expansion[step] = optic_speed
        
    def solve(self):
        self.step = 0
        for step in range(self.simlen-1):
            self.euler()
                     
    def euler (self, dt=None):
        #print self.angle_to_post[self.step], self.world_ori[self.step]
        self.step += 1
        step = self.step
        if dt is None:
            dt = self.dt
    
        phi = self.angle_to_post[step-1][0]
        ori = self.world_ori[step-1][0]
        ori_dot = 0
        
        
        s = np.array([  [ori], 
                        [ori_dot],
                     ])
        
        # angular dynamics from Reichardt and Poggio, 1976
        A = np.array(  [[0.,    1.],
                        [0.,    -self.k/self.rotmom],
                        ])
        b = np.array(  [[0],
                        [N(phi) + D(phi)/self.rotmom + r(phi)/self.rotmom],
                        ])
                        
        sdot = np.dot(A,s) + b
        new_s = s + dt*sdot  
        #print new_s[1], s[1], sdot[1]       
        self.world_ori[step] = self.world_ori[step-1][0] + dt*sdot[1]
        
        # for now, assume fly flies in the direction it is facing, with a constant velocity
        speed = np.linalg.norm(self.vel[step-1])
        
        y = np.sin(ori)
        x = np.cos(ori)
        
        print ori*180./np.pi, x, y

        vec = np.array([x,y])        
        
        new_vel = speed*vec
        new_pos = self.pos[step-1]+dt*self.vel[step-1]
        
        self.pos[step] = new_pos
        self.vel[step] = new_vel
        
        self.calc_relative_states( step )        
        #print self.angle_to_post[step], self.world_ori[step], dt*sdot[1]        
        
if __name__ == '__main__':
                
    post = Post(radius=0.1)
    fly = Fly(post, initial_pos=[-.01,-0.01], initial_vel=[0.01,0.0], initial_ori=0, dt=0.005, simlen=1000)

    fly.solve()
        
        
        
''' IDEAS

- stripe fixation response should depend on angle subtended by object, contrast, or something of the sort
- 


'''
        
