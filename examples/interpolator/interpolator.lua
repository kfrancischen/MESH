interpolator = Interpolator.new({3.0, 5.4, 5.7, 8.0},
	{{14.2, 32}, -- x, and list of y values
	 {4.6, 10},
	 {42.7, 20},
	{35.2, 40}}
	)

omega = {}
epsilon = {}
for x = 3, 8, 0.1 do
	y = interpolator:Get(x)
  table.insert(omega, x)
  table.insert(epsilon, y)
end

s = SimulationPlanar.new()
s:AddMaterial("test", omega, epsilon)