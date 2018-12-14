#! /usr/local/bin/python3.7
"""
3-D spline interpolation
(with graph drawing by matplotlib)
"""
import matplotlib.pyplot as plt
import sys
import traceback

class SplineInterpolation:
    def __init__(self, xs, ys):
        """ Initialization

        :param list xs: x-coordinate list of given points
        :param list ys: y-coordinate list of given points
        """
        self.xs, self.ys = xs, ys
        self.n = len(self.xs) - 1
        h = self.__calc_h()
        w = self.__calc_w(h)
        matrix = self.__gen_matrix(h, w)
        v = [0] + self.__gauss_jordan(matrix) + [0]
        self.b = self.__calc_b(v)
        self.a = self.__calc_a(v)
        self.d = self.__calc_d()
        self.c = self.__calc_c(v)

    def interpolate(self, t):
        """ Interpolation

        :param  float t: x-value for a interpolate target
        :return float  : computated y-value
        """
        try:
            i = self.__search_i(t)
            return self.a[i] * (t - self.xs[i]) ** 3 \
                 + self.b[i] * (t - self.xs[i]) ** 2 \
                 + self.c[i] * (t - self.xs[i]) \
                 + self.d[i]
        except Exception as e:
            raise

    def __calc_h(self):
        """ H calculation

        :return list: h-values
        """
        try:
            return [self.xs[i + 1] - self.xs[i] for i in range(self.n)]
        except Exception as e:
            raise

    def __calc_w(self, h):
        """ W calculation

        :param  list h: h-values
        :return list  : w-values
        """
        try:
            return [
                6 * ((self.ys[i + 1] - self.ys[i]) / h[i]
                   - (self.ys[i] - self.ys[i - 1]) / h[i - 1])
                for i in range(1, self.n)
            ]
        except Exception as e:
            raise

    def __gen_matrix(self, h, w):
        """ Matrix generation

        :param  list   h: h-values
        :param  list   w: w-values
        :return list mtx: generated 2-D matrix
        """
        mtx = [[0 for _ in range(self.n)] for _ in range(self.n - 1)]
        try:
            for i in range(self.n - 1):
                mtx[i][i]     = 2 * (h[i] + h[i + 1])
                mtx[i][-1]    = w[i]
                if i == 0:
                    continue
                mtx[i - 1][i] = h[i]
                mtx[i][i - 1] = h[i]
            return mtx
        except Exception as e:
            raise

    def __gauss_jordan(self, matrix):
        """ Solving of simultaneous linear equations
            with Gauss-Jordan's method

        :param  list mtx: list of 2-D matrix
        :return list   v: answers list of simultaneous linear equations
        """
        v = []
        n = self.n - 1
        try:
            for k in range(n):
                p = matrix[k][k]
                for j in range(k, n + 1):
                    matrix[k][j] /= p
                for i in range(n):
                    if i == k:
                        continue
                    d = matrix[i][k]
                    for j in range(k, n + 1):
                        matrix[i][j] -= d * matrix[k][j]
            for row in matrix:
                v.append(row[-1])
            return v
        except Exception as e:
            raise

    def __calc_a(self, v):
        """ A calculation

        :param  list v: v-values
        :return list  : a-values
        """
        try:
            return [
                (v[i + 1] - v[i])
              / (6 * (self.xs[i + 1] - self.xs[i]))
                for i in range(self.n)
            ]
        except Exception as e:
            raise

    def __calc_b(self, v):
        """ B calculation

        :param  list v: v-values
        :return list  : b-values
        """
        try:
            return [v[i] / 2.0 for i in range(self.n)]
        except Exception as e:
            raise

    def __calc_c(self, v):
        """ C calculation

        :param  list v: v-values
        :return list  : c-values
        """
        try:
            return [
                (self.ys[i + 1] - self.ys[i]) / (self.xs[i + 1] - self.xs[i]) \
              - (self.xs[i + 1] - self.xs[i]) * (2 * v[i] + v[i + 1]) / 6
                for i in range(self.n)
            ]
        except Exception as e:
            raise

    def __calc_d(self):
        """ D calculation

        :return list: c-values
        """
        try:
            return self.ys
        except Exception as e:
            raise

    def __search_i(self, t):
        """ Index searching

        :param float t: t-value
        :return  int i: index
        """
        i, j = 0, len(self.xs) - 1
        try:
            while i < j:
                k = (i + j) // 2
                if self.xs[k] < t:
                    i = k + 1
                else:
                    j = k
            if i > 0:
                i -= 1
            return i
        except Exception as e:
            raise

class Graph:
    def __init__(self, xs_0, ys_0, xs_1, ys_1):
        self.xs_0, self.ys_0, self.xs_1, self.ys_1 = xs_0, ys_0, xs_1, ys_1

    def plot(self):
        """ Graph plotting """
        try:
            plt.title("3-D Spline Interpolation")
            plt.scatter(
                self.xs_1, self.ys_1, c = "b",
                label = "interpolated points", marker = "+"
            )
            plt.scatter(
                self.xs_0, self.ys_0, c = "r",
                label = "given points"
            )
            plt.xlabel("x")
            plt.ylabel("y")
            plt.legend(loc = 2)
            plt.grid(color = "gray", linestyle = "--")
            #plt.show()
            plt.savefig("spline_interpolation.png")
        except Exception as e:
            raise


if __name__ == '__main__':
    # (N + 1) points
    X = [0.0, 2.0, 3.0, 5.0, 7.0, 8.0,18]
    Y = [0.8, 2.8, 3.2, 1.9, 4.5, 2.5,7]
    S   = 0.1        # Step for interpolation
    S_1 = 1 / S      # Inverse of S
    xs_g, ys_g = [], []  # List for graph
    try:
        # 3-D spline interpolation
        si = SplineInterpolation(X, Y)
        for x in [x / S_1 for x in range(int(X[0] / S), int(X[-1] / S) + 1)]:
            y = si.interpolate(x)
            print("{:8.4f}, {:8.4f}".format(x, y))
            xs_g.append(x)
            ys_g.append(y)
        # Graph drawing
        g = Graph(X, Y, xs_g, ys_g)
        g.plot()
    except Exception as e:
        traceback.print_exc()
        sys.exit(1)
