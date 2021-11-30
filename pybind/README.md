# A Python SpMV Script

### Run (using pybind11)

```bash
pip install pybind11

# MacOS
g++ -O3 -Wall -shared -std=c++11 -undefined dynamic_lookup `python3 -m pybind11 --includes` _csb.cpp -o _csb`python3-config --extension-suffix`

# Linux
g++ -O3 -Wall -shared -std=c++11 -fPIC `python3 -m pybind11 --includes` _csb.cpp -o _csb`python3-config --extension-suffix`

python3 test.py
```

```python
Mat = [[1, 2, 3], [0 ,4, 0], [5, 0, 0]]
Mat_CSB = _csb.CSB(Mat)
Vec = [0, 1, 2]
print(Mat_CSB.SpMV(Vec))
```

