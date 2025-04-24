from math import sqrt
from random import uniform
import typing as tp

from dichotomy import *
from copy import deepcopy
from itertools import combinations



class matrix:
    def __init__(
            self,
            n: int,
            m: int,
            values: tp.Optional[tp.List[tp.List[float]]] = None
        )-> None:
        """Основной конструктор для класса матриц (matrix)

        Args:
            n (int): кол-во строк
            m (int): кол-во стобцов
            values (tp.Optional[tp.List[tp.List[float]]]): массив значений (Optional)

        Raises:
            ValueError: в случае неккоректности n и m 
            ValueError: в случае несовпадения n и кол-ва строк в values
            ValueError: в случае несовпадения m и кол-ва столбцов в values
        """
        if (n < 1 or m < 1): raise ValueError("Incorrect matrix dimensions")

        self._det = None
        self.det_change = False
        self.n = n
        self.m = m

        if (values is None):
            self.values = [[0 for i in range(m)] for j in range(n)]
            return
        
        if (len(values) != n): raise ValueError("Incorrect values array dimension (rows)")

        for i in range(len(values)):
            if len(values[i]) != m: raise ValueError("Incorrect values array dimension (columns)")

        self.values = deepcopy(values)
    
    def __getitem__(self, index: tuple[int, int]) -> float:
        """Получить элемент матрицы по координатам

        Args:
            index (tuple[int, int]): координаты

        Raises:
            ValueError: в случае неккоректных координат

        Returns:
            float: элемент по координатам 
        """
        if index[0] < 0 or index[0] >= self.n or index[1] < 0 or index[1] >= self.m:
            raise ValueError("Incorrect coordinates")
        
        return self.values[index[0]][index[1]]

    def __setitem__(self, index: tuple[int, int], value: float) -> None:
        """Присвоить элемент матрице по координатам

        Args:
            index (tuple[int, int]): координаты
            value (float): значение

        Raises:
            ValueError: в случае неккоректных координат
        """
        self.det_change = True
        if index[0] < 0 or index[0] >= self.n or index[1] < 0 or index[1] >= self.m:
            raise ValueError("Incorrect coordinates")
        
        self.values[index[0]][index[1]] = value

    def __add__(self, other) -> tp.Self:
        """Сложить две матрицы соответствующих типов или единичную матрицу 
        соответствующего размера, заданную скаляром

        Args:
            other (_type_): матрица или скаляр

        Raises:
            ValueError: в случае несовпадения размерностей матриц
            ValueError: в случае неподходящего типа other

        Returns:
            tp.Self: результирующая матрица (sum matrix)
        """
        res = matrix(self.n, self.m, self.values)

        if isinstance(other, (int, float)):
            for i in range(min(self.m, self.n)): res[i, i] += other

        elif isinstance(other, matrix):
            if self.n != other.n or self.m != other.m: raise ValueError("Incompatible matrix dimensions")

            for i in range(self.n):
                for j in range(self.m):
                    res[i, j] += other[i, j]
        else: raise ValueError("Incompatible types")

        return res

    def __radd__(self, other) -> tp.Self:
        return self.__add__(other)
    
    def __sub__(self, other) -> tp.Self:
        return self.__add__((-1) * other)

    def __mul__(self, other) -> tp.Self:
        """Умножение матрицы на скаляр и другую матрицу

        Args:
            other (_type_): скаляр или matrix

        Raises:
            ValueError: в случае несовместимых размерностей
            ValueError: в случае несовместимости типов

        Returns:
            tp.Self: результирующая матрица (mult matrix)
        """
        if isinstance(other, (int, float)):
            res = matrix(self.n, self.m, self.values)
            for i in range(self.n):
                for j in range(self.m):
                    res[i, j] *= other
        elif isinstance(other, matrix):
            if self.m != other.n: raise ValueError("Incompatible matrix dimensions")
            res = matrix(self.n, other.m)
            for i in range(self.n):
                for j in range(other.m):
                    s = 0
                    for k in range(self.m): s += self[i, k] * other[k, j]
                    res[i, j] = s
        else: raise ValueError("Incompatible types")

        return res
    
    def __rmul__(self, other):
        return self.__mul__(other)
    
    def __str__(self) -> str:
        mtrx = ''
        for i in range(self.n):
            for j in range(self.m):
                mtrx += f'{round(self[i, j], 4)} '
            mtrx += '\n'
        
        return mtrx

    def concat(self, other: tp.Self, axis: tp.Literal[0, 1]) -> tp.Self:
        if not axis in [0, 1]: raise ValueError("Incorrect axis chosen")
        if (axis == 0):
            if (self.m != other.m): raise ValueError("Dimensions (cols) do not correlate!")
            return self.T.concat(other=other.T, axis=1).T
        else:
            if (self.n != other.n): raise ValueError("Dimensions (rows) do not correlate!")
            res = matrix(self.n, self.m + other.m, values=[self.values[i] + other.values[i] for i in range(self.n)])
            return res
            

    def create_minor(self, I: tuple, J: tuple) -> tp.Self:
        """Создать минор, вычеркивая определенные строки и столбцы

        Args:
            I (tuple): номера вычеркиваемых строк
            J (tuple): номера вычеркиваемых столбцов

        Returns:
            tp.Self: матрица, соответствующая минору
        """
        res = matrix(self.n - len(I), self.m - len(J))

        res.values = [[self[i, j] for j in range(self.m) if not j in J] for i in range(self.n) if not i in I]

        return res

    def calc_det(self) -> float:
        """Вычислить детерминант матрицы

        Raises:
            ValueError: В случае неквадратной матрицы

        Returns:
            float: детерминант
        """
        if self.n != self.m:
            raise ValueError("The matrix is not square")
        matrix = deepcopy(self.values)
        det = 1.0

        for i in range(self.n):
            # Поиск максимального элемента в текущем столбце для выбора ведущего элемента
            max_row = i
            for k in range(i, self.n):
                if abs(matrix[k][i]) > abs(matrix[max_row][i]):
                    max_row = k

            # Если ведущий элемент ноль, то определитель равен нулю
            if matrix[max_row][i] == 0:
                return 0

            # Меняем строки местами, если нужно
            if max_row != i:
                matrix[i], matrix[max_row] = matrix[max_row], matrix[i]
                det *= -1  # Меняем знак определителя из-за перестановки строк

            # Прямой ход метода Гаусса
            for j in range(i + 1, self.n):
                factor = matrix[j][i] / matrix[i][i]
                for k in range(i, self.n):
                    matrix[j][k] -= factor * matrix[i][k]

            # Умножаем определитель на диагональный элемент
            det *= matrix[i][i]
        return det
        

    @property
    def det(self) -> float:
        if (self.n != self.m): raise AttributeError("Can't compute det for non-square matrix")
        if (self._det is None or self.det_change):
            self._det = self.calc_det()
            self.det_change = False
        return self._det    

    @property
    def T(self) -> tp.Self:
        """Вернуть транспонированную матрицу

        Returns:
            tp.Self: транспонированная матрица
        """
        res = matrix(self.m, self.n)
        for i in range(self.m):
            for j in range(self.n):
                res[i, j] = self[j, i]
        
        return res
    
    @property
    def symmetry(self) -> bool:
        """Проверить матрицу на симметричность

        Returns:
            bool: флаг симметричности
        """
        if (self.n != self.m): return False

        for i in range(self.n):
            for j in range(i + 1, self.m):
                if self[i, j] != self[j, i]: return False

        return True 

    @property
    def l2_norm(self) -> float:
        return sqrt(sum(self.values[i][j] ** 2 for i in range(self.n) for j in range(self.m)))

class polynomial:
    def __init__(
        self,
        coefs: tp.List[float]
    ) -> None:
        """Конструктор для класса полиномов

        Args:
            coefs (tp.List[float]): коэффициенты полинома (начиная с высшего)
        """

        self.n = len(coefs) - 1
        self.coefs = coefs.copy()

    def __getitem__(self, coord: int) -> float:
        if (coord < 0): raise ValueError("Invalid coord")
        if (coord > self.n): return 0.0

        return self.coefs[coord]
    
    def __setitem__(self, coord: int, value: float) -> None:
        if (coord < 0): raise ValueError("Invalid coord")

        if (coord <= self.n): self.coefs[coord] = value
        else:
            coefs_new = [0 for _ in range(coord + 1)]
            for i in range(self.n + 1): coefs_new[i] += self.coefs[i]

            coefs_new[coord] = value
            self.coefs = coefs_new
            self.n = coord

    def value(self, x: float) -> float:
        """Вычислисть значение полинома в точке

        Args:
            x (float): точка

        Returns:
            float: значение полинома
        """
        x_c = 1
        y = 0
        for i in range(self.n + 1):
            y += x_c * self.coefs[i]
            x_c *= x
        
        return y
    
    def __str__(self) -> str:
        pol = f't^{self.n}'
        for i in range(self.n - 1, 0, -1):
            c = self.coefs[i]
            if c != 0: pol += ((' + ' if c > 0 else ' - ') + f'{round(abs(c), 3)}t^{i}')
        c = self.coefs[0]
        if c != 0: pol += ((' + ' if c > 0 else ' - ') + f'{round(abs(c), 3)}')
        
        return pol

def find_minor_sum(A: matrix, k: int) -> float:
    """Найти сумму главных миноров k-ого порядка

    Args:
        A (matrix): матрица
        k (int): порядок

    Raises:
        ValueError: в случае неквадратной матрицы
        ValueError: в случае некорректного порядка

    Returns:
        float: сумма главных миноров
    """
    if (A.n != A.m): raise ValueError("Non-square matrix passed")
    if (k < 0 or k > A.n): raise ValueError("Invalid k value")

    if (k == A.n): return 1

    sum_minors = 0
    for c in combinations(range(A.n), k): sum_minors += A.create_minor(c, c).det

    return sum_minors



def calc_char_pol(C: matrix) -> polynomial:
    """Вычислить коэффициенты характеристического многочлена

    Args:
        C (matrix): матрица

    Raises:
        ValueError: в случае неквадратной матрицы

    Returns:
        polynomial: характеристический многочлен
    """
    if (C.n != C.m): raise ValueError("Non-square matrix passed")

    coefs = [0 for _ in range(C.n + 1)]

    for i in range(C.n + 1):
        coefs[i] = find_minor_sum(C, i) * ((-1) ** ((C.n - i) % 2))
    
    return polynomial(coefs)


# EASY

def gauss_solver(A: matrix, b: matrix, zero_err: float = 10e-4) -> tuple[matrix, tp.List[matrix]]:
    """Решить СЛАУ методом Гаусса

    Args:
        A (matrix): матрица коэффициентов
        b (matrix): вектор свободных членов
        zero_err (float, optional): порог нуля. Defaults to 10e-4.

    Raises:
        NotImplementedError: так как в данной лабораторной работе не требуется 
        находить решение неквадратной СЛАУ, данный метод специализирован под
        квадратные 
        ValueError: в случае несоответствия размерностей матрицы коэффициентов и свободных членов
        ValueError: в случае, если вектор свободных членов имеет неподходящую размерность
        ValueError: в случае несовместности системы

    Returns:
        tuple[matrix, tp.List[matrix]]: кортеж частного решения и ФСР однородной системы
    """
    if (A.n != A.m): raise NotImplementedError("Gauss solver not yet implemented for non-square matrices")

    if (A.n != b.n): raise ValueError("Coefs and free terms dimensions do not correlate!")
    if (b.m != 1): raise ValueError("Free term vector is dimensity 1 on axis 1")

    ext = A.concat(b, axis=1)
    n = A.n
    for i in range(n):
        max_row = i
        for j in range(i, n):
            if abs(ext[j, i]) > abs(ext[max_row, i]): max_row = j
        
        if ext[max_row, i] == 0: continue

        ext.values[max_row], ext.values[i] = ext.values[i], ext.values[max_row]

        p = ext[i, i]

        for j in range(i + 1, n):
            factor = ext[j, i] / p
            for k in range(i, n + 1):
                ext[j, k] -= ext[i, k] * factor
    
    X = [0] * n
    basis_vec = []
    free_vars = []

    for i in range(n - 1, -1, -1):
        r = ext[i, n] - sum(ext[i, j] * X[j] for j in range(i + 1, n))
        if abs(ext[i, i]) < zero_err:
            if abs(r) > zero_err:
                raise ValueError("The system is inconsistent")
            free_vars.append(i)
            continue
        X[i] = r / ext[i, i]
    X = matrix(n, 1, [[x] for x in X])

    for var in free_vars:
        vec = [0] * n
        vec[var] = 1
        for i in range(n - 1, -1, -1):
            if i in free_vars: continue
            vec[i] -= sum(ext[i, j] * vec[j] for j in range(i + 1, n)) / ext[i, i]
        basis_vec.append(matrix(n, 1, [[x] for x in vec]))
    
    return X, basis_vec

def center_data(X: matrix) -> matrix:
    """Отцентрировать данные по среднему в столбце

    Args:
        X (matrix): матрица данных

    Returns:
        matrix: отцентрированная матрица
    """
    X_c = matrix(X.n, X.m, X.values)
    means = [0 for i in range(X.m)]
    for j in range(X.m):
        s = 0.0
        for i in range(X.n):
            s += X[i, j]
        s /= X.n
        means[j] = s
    for i in range(X.n):
        for j in range(X.m):
            X_c[i, j] -= means[j]
    
    return X_c

def covariance_matrix(X_c: matrix) -> matrix:
    """Вычислить матрицу ковариаций

    Args:
        X_c (matrix): отцентрированная матрица данных

    Returns:
        matrix: матрица ковариаций
    """
    return X_c.T*X_c*(1/(X_c.n - 1))

# NORMAL

def find_eigenvalues(C: matrix, err: float = 1e-6, enable_gershgorin: bool = False) -> tp.List[float]:
    """Найти собственные значения

    Args:
        C (matrix): квадратная матрица
        err (float, optional): ошибка вычисления. Defaults to 1e-6.
        enable_gershgorin (bool, optional): использовать ли теорему Гершгорина для 
        нахождения промежутков, содержащих значения. Defaults to False.

    Returns:
        tp.List[float]: массив собственных значений
    """
    if (C.n != C.m): raise ValueError("Passed down a non-square matrix")

    char_pol = calc_char_pol(C)

    if enable_gershgorin:
        roots = []
        for i in range(C.n):
            r = 0
            for j in range(C.m):
                r += abs(C[i, j]) if j != i else 0
            roots += find_root_advanced(C[i, i] - r, C[i, i] + r, char_pol.value, err=err)
    else:
        roots = find_root_advanced(-100, 100, char_pol.value, err=err)
    
    res = []
    for i in roots:
        present = False
        for j in res:
            if abs(i - j) <= 2 * err:
                present = True
                break
        if not present:
            res.append(i)
    res = sorted(res)
    res.reverse()
    return res

def find_eigenvectors(C: matrix, eigenvalues: tp.List[float], zero_err: None | float = None) -> tuple[tp.List[matrix], tp.List[float]]:
    """Найти собственные векторы матрицы

    Args:
        C (matrix): матрица
        eigenvalues (tp.List[float]): массив собственных значений (спектр)

    Returns:
        tuple[tp.List[matrix], tp.List[float]] массив собственных векторов и обновленный (с учтенными кратностями) спектр
    """
    if zero_err is None: zero_err = 1e-4
    eigenvectors = []
    spectrum = []
    b = matrix(C.n, 1, [[0] for _ in range(C.n)])
    for lam in eigenvalues:
        for vector in gauss_solver(C - lam, b, zero_err)[1]:
            vector = vector * (1 / vector.l2_norm)
            eigenvectors.append(vector)
            spectrum.append(lam)

    return eigenvectors, spectrum

def explained_variance_ratio(eigenvalues: tp.List[float], k: int) -> float:
    """Вычислить долю объясненной дисперсии

    Args:
        eigenvalues (tp.List[float]): собственные значения
        k (int): число компонент

    Returns:
        float: доля объясненной дисперсии
    """
    return sum(eigenvalues[:k])/sum(eigenvalues)

# HARD

def pca(X: matrix, k: int | None = None, auto_select_k_threshold: tp.Optional[float] = None,
        use_zero_vals: bool = False, zero_err: None | float = None, enable_gershgorin: bool = False) -> tuple[matrix, float]:
    """Алгоритм PCA

    Args:
        X (matrix): матрица данных
        k (None | int, optional): кол-во компонент
        auto_select_k_threshold (None | float): порог дисперсии (присвоить float для автом. выбора k)
        zero_vals (bool, optional): вставлять ли в результирующие данные нулевые столбцы
        zero_err (float, optional): порог нуля (для find_eigenvectors)
        enable_gershgorin (bool, optional): использовать ли теорему Гершгорина для нахождения собственных знач.

    Raises:
        ValueError: в случае невалидности выбранного k
        ValueError: в случае, если в оба параметра k и auto_select_k_threshold переданы None
        ArithmeticError: в случае ошибки вычисления собственных векторов (change find_eigenvectors err)

    Returns:
        tuple[matrix, float]: X_proj и дисперсия
    """
    if not k is None and (k < 1 or k > X.m): raise ValueError("Invalid k (must be in the 1 to X.m range)")

    if k is None and auto_select_k_threshold is None: raise ValueError("No instruction on k given")

    X_c = center_data(X)
    cov = covariance_matrix(X_c)

    eigenvals = find_eigenvalues(cov, enable_gershgorin=enable_gershgorin)
    eigenvecs, spectrum = find_eigenvectors(cov, eigenvals, zero_err)

    if not auto_select_k_threshold is None:
        k = auto_select_k(spectrum, threshold=auto_select_k_threshold)
    if (len(eigenvecs) != X.m): raise ArithmeticError("Can't compute eigenvecs!")

    Vk = eigenvecs[0]
    for i in range(1, k):
        Vk = Vk.concat(eigenvecs[i], axis=1)
    if use_zero_vals:
        zero_col = matrix(X.m, 1)
        for i in range(k, X.m):
            Vk = Vk.concat(zero_col, axis=1)

    X_proj = X_c * Vk

    return X_proj, explained_variance_ratio(spectrum, k)

def reconstruction_error(X_or: matrix, X_rec: matrix) -> float:
    """Вычислить меру ошибки реконструкции данных

    Args:
        X_or (matrix): оригинальная матрица данных
        X_rec (matrix): реконструированная

    Returns:
        float: MSE
    """
    err = 0.0
    for i in range(X_or.n):
        for j in range(X_or.m):
            err += (X_or[i, j] - X_rec[i, j]) ** 2
    return err / (X_or.n * X_or.m)

def auto_select_k(spectrum: tp.List[float], threshold: float = 0.95) -> int:
    """Выбрать кол-во компонент на основе порога дисперсии

    Args:
        spectrum (tp.List[float]): спектр матрицы (отсортирован по убыванию)
        threshold (float, optional): порог. Defaults to 0.95.

    Returns:
        int: кол-во компонент
    """
    sum_lam_k = 0
    k = 0
    sum_lam = sum(spectrum)
    while (sum_lam_k/sum_lam < threshold and k < len(spectrum)):
        sum_lam_k += spectrum[k]
        k += 1
    
    return k

def reconstruct_data(X: matrix, X_proj: matrix, Vk: matrix) -> matrix:
    means = [0] * X.m
    for i in range(X.n):
        for j in range(X.m):
            means[j] += X[i, j]
    for j in range(X.m): means[j] /= X.n

    X_rec = X_proj * Vk.T

    for i in range(X.n):
        for j in range(X.m):
            X_rec[i, j] += means[j]
    
    return X_rec

def handle_missing_values(X: matrix) -> matrix:
    """Обработать пропущенные значения (заменить на среднее признака)

    Args:
        X (matrix): матрица данных

    Returns:
        matrix: 
    """
    means = [0] * X.m
    X_fix = matrix(X.n, X.m, X.values)
    for j in range(X.m):
        not_nan = 0
        for i in range(X.n):
            if X[i, j] != X[i, j]: continue 
            not_nan += 1
            means[j] += X[i, j]
        means[j] /= not_nan
    
    for i in range(X.n):
        for j in range(X.m):
            if X_fix[i, j] != X_fix[i, j]: X_fix[i, j] = means[j]
    
    return X_fix

def add_noise_and_compare(X: matrix, noise_level: float = 0.1, k: int | None = None,
                            auto_select_k_threshold: float | None = None) -> tuple[tuple[int, float, float], tuple[int , float, float]]:
    """Сравнить метрики после добавления шума в датасет

    Args:
        X (matrix): матрица данных
        noise_level (float, optional): уровень шума. Defaults to 0.1.
        k (int | None, optional): кол-во компонент. Defaults to None.
        auto_select_k_threshold (float | None, optional): порог дисперсии для автом. выбора k. Defaults to None.

    Raises:
        ValueError: в случае отрицательного noise level
        ValueError: в случае, если в оба параметра k и auto_select_k_threshold переданы None

    Returns:
        tuple[tuple[int, float, float], tuple[int , float, float]]: 
        ((кол-во компонент изнач. матрицы, дисперсия, MSE), (те же метрики, но для матрицы с шумом)
    """
    if (noise_level < 0): raise ValueError("Incoherent noise level")
    if k is None and auto_select_k_threshold is None: raise ValueError("No instructions on k given")

    sigma = [0] * X.m
    X_c = center_data(X)
    for i in range(X.n):
        for j in range(X.m):
            sigma[j] += X_c[i, j] ** 2
    for j in range(X.m): sigma[j] = sqrt(sigma[j])

    X_n = matrix(X.n, X.m, values=X.values)
    for i in range(X.n):
        for j in range(X.m):
            X_n[i, j] += uniform(-noise_level * sigma[j], noise_level * sigma[j])
    
    cov = covariance_matrix(X_c)
    eigenvals = find_eigenvalues(cov, enable_gershgorin=True)
    eigenvecs, spectrum = find_eigenvectors(cov, eigenvals)

    X_n_c = center_data(X_n)
    cov_n = covariance_matrix(X_n_c)
    eigenvals_n = find_eigenvalues(cov_n, enable_gershgorin=True)
    eigenvecs_n, spectrum_n = find_eigenvectors(cov_n, eigenvals_n, zero_err=1e-3)

    if not auto_select_k_threshold is None:
        k = auto_select_k(spectrum, auto_select_k_threshold)
        k_n = auto_select_k(spectrum_n, auto_select_k_threshold)
    else:
        k_n = k
    
    X_proj, disp = pca(X, k, use_zero_vals=True, enable_gershgorin=True)
    X_n_proj, disp_n = pca(X_n, k_n, use_zero_vals=True, zero_err=1e-3, enable_gershgorin=True)

    Vk = eigenvecs[0]
    for i in range(1, k):
        Vk = Vk.concat(eigenvecs[i], axis=1)
    zero_col = matrix(X.m, 1)
    for i in range(k, X.m):
        Vk = Vk.concat(zero_col, axis=1)

    Vk_n = eigenvecs_n[0]
    for i in range(1, k):
        Vk_n = Vk_n.concat(eigenvecs_n[i], axis=1)
    for i in range(k, X.m):
        Vk_n = Vk_n.concat(zero_col, axis=1)
    
    X_rec = reconstruct_data(X, X_proj, Vk)
    X_n_rec = reconstruct_data(X_n, X_n_proj, Vk_n)

    rec_err = reconstruction_error(X, X_rec)
    rec_n_err = reconstruction_error(X_n, X_n_rec)

    return ((k, disp, rec_err), (k_n, disp_n, rec_n_err))