
import numpy as np


class A:

    def __init__(self, data):
        self.data = data

    @classmethod
    def from_arg(cls, arg):
        data = Test._private_func(arg)
        print(cls._normalize(np.array([1, 2, 1])))
        return cls(data)
    @classmethod
    def _private_func(cls, arg):
        return arg + 1

    @classmethod
    def _normalize(cls, vector):
        """returns a normalized vector. Input, output is an numpy array."""
        norm = np.linalg.norm(vector)
        if norm == 0:
            return vector
        return vector / norm

    def add_to_b(self):
        a_list = [1, 1, 1]
        self.data[1].add_numbers(a_list)

class B:
    def __init__(self, data):
        self.data = data
        self.numbers = []

    def add_numbers(a_list):
        for num in a_list:
            self.numbers.append(num)



def main():
    #test1 = Test(1)
    #print(test1.data)
    #test2 = Test.from_arg(0)
    #print(test2.data)

    x = B(1)
    y = B(2)
    print (y.numbers
    a = A([x, y])
    a.add_to_b()
    print(y.numbers


main()
