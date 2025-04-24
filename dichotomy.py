import typing as tp

def dichotomy_search(a: float, b: float, f: tp.Callable, err: float) -> float:
    epsilon = err / 2
    while (b - a) / 2 > err:
        y_k = (a + b - epsilon) / 2
        z_k = (a + b + epsilon) / 2
        if f(y_k) <= f(z_k):
            b = z_k
        else:
            a = y_k
    return (a + b) / 2

def find_root(a: float, b: float, f: tp.Callable, err: float) -> tp.Optional[float]:
    if f(a)*f(b) > 0:
        return None  
    while (b - a) > err:
        c = (b + a) / 2
        if f(a)*f(c) > 0:
            a = c
        else:
            b = c
    return (b + a) / 2

def find_root_advanced(a: float, b: float, f: tp.Callable, err: float,
                       part_range: tp.Optional[float] = None, err_dichotomy: tp.Optional[float] = None) -> tp.List[float]:
    
    def mod_f(x: float):
        return abs(f(x))

    if part_range is None: part_range = err * 100
    roots = []
    a_1 = a
    b_1 = a + part_range
    while b_1 <= b:
        if f(a_1)*f(b_1) <= 0:
            roots.append(find_root(a_1, b_1, f, err))

        elif not err_dichotomy is None:
            prediction = dichotomy_search(a=a_1, b=b_1, f=mod_f, err=err)

            if mod_f(prediction) < err_dichotomy:
                roots.append(prediction)

        a_1 = b_1
        b_1 = b_1 + part_range
    return roots