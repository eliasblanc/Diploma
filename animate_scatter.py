import matplotlib.pyplot as plt
import matplotlib.animation as animation

f = open('royal.txt', 'r')
tdata = f.readlines()

plt.rcParams["figure.figsize"] = [5, 5]
plt.rcParams["figure.autolayout"] = True

data = []
sol = []

for line in tdata:
	points = line.split()
	data.append(points)
	sol.append('blue')

data = [list(map(float, sublist)) for sublist in data]

fig, ax = plt.subplots()
marker_size = 0.3

index = 0

def animate(i):
	fig.clear()
	ax = fig.add_subplot(111, aspect='equal', autoscale_on=False, xlim=(-1, 1), ylim=(-1, 1))
	ax.set_xlim(-1, 1)
	ax.set_ylim(-1, 1)
	global index
	a = i // 10
	n = (a + 1) * 9
	for j in range(n):
		s = ax.scatter(data[index][0], data[index][1], s=marker_size, c=sol[index])#, cmap="RdBu_r", marker="o", edgecolor='black')
		index += 1

plt.grid(b=None)
ani = animation.FuncAnimation(fig, animate, interval=30, frames=100)

ani.save('animation.gif', writer='pillow')
