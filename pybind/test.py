import _csb

def SpMV(Matrix, Vector):
	CSB = _csb.CSB(Matrix)
	return CSB.SpMV(Vector)

# a = [[1, 2, 3], [0, 4, 0], [5, 0, 0]]
# b = [1, 2, 3]
# print(SpMV(a, b))

n, m = input().split()
mat = []
for i in range(int(n)):
	mat.append([float(i) for i in input().split()])
vec = [float(i) for i in input().split()]
print(SpMV(mat, vec))