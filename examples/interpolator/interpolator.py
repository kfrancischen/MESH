from MESH import Interpolator, SimulationPlanar

interpolator = Interpolator( ((3.0, (14.2, 32)), (5.4, (4.6, 10)), (5.7, (42.7, 20)),(8.0, (35.2, 40))) )
data = ()
for i in range(30, 81, 1):
  eps = interpolator.Get(i * 0.1)
  data += ((i*0.1, eps),)

print data
s = SimulationPlanar()
s.AddMaterial( "test", data )