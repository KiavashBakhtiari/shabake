# %% [markdown]
# شبیه‌ساز کامل BCH (باینری) — Jupyter-compatible Python script
#
# این نوت‌بوکِ کد (فرمت .py با سلول‌های Jupyter) یک شبیه‌ساز آموزشی و قابل اجرا
# برای کدهای BCH (narrow-sense binary BCH) فراهم می‌کند. قابلیت‌ها:
# - ساخت میدان GF(2^m) با جداول log/antilog
# - محاسبهٔ چندجمله‌ای‌های مینیمال با استفاده از کلاس‌های کونژوگاسیون
# - ساخت چندجمله‌ای مولد g(x)
# - انکُدینگ سیستماتیک
# - دیکُدینگ: محاسبهٔ سیندروم‌ها، الگوریتم Berlekamp–Massey، Chien search و اصلاح خطا
# - شبیه‌سازی کانال BSC و ارزیابی عملکرد (FER/BER)
#
# توضیح: این پیاده‌سازی برای اهداف آموزشی طراحی شده و برای m تا حدوداً 6 یا 7
# (n = 2^m - 1 تا 127) مناسب است. برای موارد بسیار بزرگ بهینه‌سازی لازم است.

# %%
# وابستگی‌ها
import numpy as np
import random
from typing import List, Tuple

# %% [markdown]
# بخش ۱ — پیاده‌سازی GF(2^m) ساده (نمایش عناصر به صورت صحیح‌های 0..2^m-1)

# %%
class GF2m:
    def __init__(self, m: int, prim_poly: int):
        """
        m: درجه میدان
        prim_poly: چندجمله‌ای مقدم بصورت باینری (مثلاً 0b10011 برای x^4 + x + 1)
        """
        assert m >= 1
        self.m = m
        self.prim = prim_poly
        self.n = (1 << m) - 1
        # جداول:
        self.exp = [0] * (2 * self.n)
        self.log = [-1] * (self.n + 1)
        self._build_tables()

    def _build_tables(self):
        x = 1
        for i in range(self.n):
            self.exp[i] = x
            self.log[x] = i
            x <<= 1
            if x & (1 << self.m):
                x ^= self.prim
        for i in range(self.n, 2 * self.n):
            self.exp[i] = self.exp[i - self.n]

    def add(self, a: int, b: int) -> int:
        return a ^ b

    def mul(self, a: int, b: int) -> int:
        if a == 0 or b == 0:
            return 0
        return self.exp[self.log[a] + self.log[b]]

    def inv(self, a: int) -> int:
        if a == 0:
            raise ZeroDivisionError("inverse of 0")
        return self.exp[self.n - self.log[a]]

    def pow(self, a: int, e: int) -> int:
        if a == 0:
            return 0
        return self.exp[(self.log[a] * e) % self.n]

    def element_to_bits(self, a: int) -> List[int]:
        return [(a >> (self.m - 1 - i)) & 1 for i in range(self.m)]

    def bits_to_element(self, bits: List[int]) -> int:
        v = 0
        for b in bits:
            v = (v << 1) | (b & 1)
        return v

# %% [markdown]
# تابع‌های کمکی برای چندجمله‌ای‌ها (دودویی)

# %%

def trim_poly(p: List[int]) -> List[int]:
    while len(p) > 1 and p[0] == 0:
        p.pop(0)
    return p


def poly_mul_bin(a: List[int], b: List[int]) -> List[int]:
    res = [0] * (len(a) + len(b) - 1)
    for i, ai in enumerate(a):
        if ai:
            for j, bj in enumerate(b):
                res[i + j] ^= bj
    return trim_poly(res)


def poly_add_bin(a: List[int], b: List[int]) -> List[int]:
    # طول‌ها را برابر کن
    la = len(a); lb = len(b)
    if la < lb:
        a = [0] * (lb - la) + a
    if lb < la:
        b = [0] * (la - lb) + b
    return trim_poly([ (ai ^ bi) for ai, bi in zip(a, b) ])


def poly_divmod_bin(dividend: List[int], divisor: List[int]) -> Tuple[List[int], List[int]]:
    # برگرداندن (quotient, remainder) در GF(2)
    A = dividend.copy()
    A = trim_poly(A)
    D = trim_poly(divisor)
    if D == [0]:
        raise ZeroDivisionError("polynomial division by zero")
    n = len(D)
    q = [0] * max(len(A) - n + 1, 1)
    while len(A) >= n:
        if A[0] == 1:
            qpos = len(A) - n
            q[qpos] = 1
            for i in range(n):
                A[i] ^= D[i]
        A.pop(0)
        A = trim_poly(A)
        if not A:
            A = [0]
    return trim_poly(q), trim_poly(A)

# %% [markdown]
# بخش ۲ — چندجمله‌ای‌های مینیمال و چندجمله‌ای مولد
# روش: برای ریشه alpha^i، کلاس کونژوگاسیون را محاسبه کن؛
# سپس چندجمله‌ای مینیمال را با ضرب عوامل (x + alpha^{e}) در فضای میدان بساز.
# سپس ضرایب میدان را باید عناصر GF(2) (یعنی 0 یا 1) باشند — در حالت narrow-sense این اتفاق می‌افتد.

# %%

def minimal_polynomial(alpha_pow: int, gf: GF2m) -> List[int]:
    """
    برمی‌گرداند ضرایب چندجمله‌ای مینیمال به صورت لیست دودویی (MSB..LSB)
    برای ریشه alpha^{alpha_pow}.
    """
    n = gf.n
    # کلاس کونژوگاسیون
    conj = []
    seen = set()
    x = alpha_pow % n
    while x not in seen:
        seen.add(x)
        conj.append(x)
        x = (x * 2) % n
    # ساخت چندجمله‌ای مینیمال به‌صورت ضرب عوامل (x + alpha^{e}) با ضرایب در GF(2^m)
    # نگهداری چندجمله‌ای به‌صورت لیست ضرایب (اولین خانه ضریب x^{deg}) اما ضرایب به‌عنوان عنصر میدان نگهداری می‌شوند
    poly_field = [1]  # مقدار 1 (coeff for x^0 initially)
    # We'll build polynomial in ascending order (const..highest) for ease, then reverse later
    poly_field = [1]
    for e in conj:
        root = gf.exp[e]  # alpha^e (field element)
        # factor is (x + root) => coefficients [root (const), 1 (x)] in ascending order
        # multiply poly_field by factor in field-coefficient domain
        new = [0] * (len(poly_field) + 1)
        for i, coeff in enumerate(poly_field):
            # multiply coeff * root goes to position i (const term)
            # new[i] += coeff * root
            new[i] ^= gf.mul(coeff, root) if coeff != 0 and root != 0 else (coeff & root)
            # multiply coeff * x goes to position i+1 (coeff stays as coeff)
            new[i + 1] ^= coeff
        poly_field = new
    # poly_field currently ascending order (const..x^L)
    # convert to descending order and ensure coefficients are in GF(2)
    poly_field_desc = poly_field[::-1]
    # map field coefficients to bits (they should be 0 or 1 for minimal polynomial over GF(2))
    bin_poly = []
    for a in poly_field_desc:
        if a == 0:
            bin_poly.append(0)
        elif a == 1:
            bin_poly.append(1)
        else:
            # Ideally should not happen; however due to representation it might produce a field element
            # that nonetheless belongs to GF(2) subfield. We check if its binary representation is 1.
            if a == 1:
                bin_poly.append(1)
            else:
                # Fallback: try to see if coefficient expressed in GF(2) basis has only lower bit 1 and others 0
                # Convert to bit vector and check if it's equal to 1
                bits = gf.element_to_bits(a)
                if bits == [0] * (gf.m - 1) + [1] or bits == [1]:
                    bin_poly.append(1)
                else:
                    # اگر به‌طور جدی ضریب در GF(2) نیست: اعلام خطا (برای ایمنی)، ولی برای اکثر پارامترها نباید رخ دهد
                    raise ValueError(f"Coefficient {a} of minimal polynomial is not in GF(2)")
    return trim_poly(bin_poly)


def bch_generator_poly(m: int, t: int, prim_poly: int) -> List[int]:
    gf = GF2m(m, prim_poly)
    n = gf.n
    polys = []
    used = set()
    for i in range(1, 2 * t + 1):
        if i in used:
            continue
        # compute conjugacy class for i
        cls = []
        x = i % n
        while x not in cls:
            cls.append(x)
            x = (x * 2) % n
        for e in cls:
            used.add(e)
        # minimal polynomial for this class
        p = minimal_polynomial(i, gf)
        polys.append(p)
    # multiply all binary minimal polynomials
    g = [1]
    for p in polys:
        g = poly_mul_bin(g, p)
    return trim_poly(g)

# %% [markdown]
# بخش ۳ — انکُدینگ سیستماتیک

# %%

def systematic_bch_encode(msg_bits: List[int], g: List[int]) -> List[int]:
    k = len(msg_bits)
    deg_g = len(g) - 1
    # message polynomial is MSB..LSB; append deg_g zeros
    padded = msg_bits + [0] * deg_g
    _, remainder = poly_divmod_bin(padded, g)
    # ensure remainder length == deg_g
    rem = remainder
    if len(rem) < deg_g:
        rem = [0] * (deg_g - len(rem)) + rem
    codeword = msg_bits + rem
    return codeword

# %% [markdown]
# بخش ۴ — محاسبهٔ سیندروم‌ها (سیندروم j برابر است با c(alpha^{j}))

# %%

def compute_syndromes(codeword: List[int], t: int, gf: GF2m) -> List[int]:
    n = gf.n
    # codeword assumed MSB->LSB (coeff x^{n-1} .. x^0)
    synd = []
    for j in range(1, 2 * t + 1):
        a = gf.exp[j]
        s = 0
        for coef in codeword:
            s = gf.mul(s, a) ^ coef
        synd.append(s)
    return synd

# %% [markdown]
# بخش ۵ — الگوریتم Berlekamp–Massey روی میدان GF(2^m)

# %%

def berlekamp_massey(synd: List[int], gf: GF2m) -> List[int]:
    # Returns sigma polynomial coefficients C[0..L] where sigma(x) = C[0] + C[1] x + ... with C[0]=1
    N = len(synd)
    C = [1] + [0] * N
    B = [1] + [0] * N
    L = 0
    m = 1
    b = 1
    for n in range(N):
        # discrepancy d = s[n] + sum_{i=1..L} C[i] * s[n-i]
        d = synd[n]
        for i in range(1, L + 1):
            if C[i] != 0 and (n - i) >= 0:
                d = d ^ gf.mul(C[i], synd[n - i])
        if d == 0:
            m += 1
        else:
            T = C.copy()
            coef = gf.mul(d, gf.inv(b))
            # C = C - coef * x^m * B  (subtraction = addition in GF(2^m) when coefficients are field elements -> xor)
            for i in range(len(B)):
                if B[i] != 0:
                    C[i + m] ^= gf.mul(coef, B[i])
            if 2 * L <= n:
                L_new = n + 1 - L
                B = T
                b = d
                L = L_new
                m = 1
            else:
                m += 1
    # trim to length L
    C = C[:L + 1]
    return C

# %% [markdown]
# بخش ۶ — Chien search برای یافتن ریشه‌های sigma(x) و اصلاح بیت‌ها (برای کد دودویی، مقدار خطا همیشه 1)

# %%

def chien_search_and_correct(codeword: List[int], sigma: List[int], gf: GF2m) -> Tuple[List[int], List[int]]:
    n = gf.n
    L = len(sigma) - 1
    error_positions = []
    cw = codeword.copy()
    # sigma is C[0] + C[1] x + ... + C[L] x^L with C[0]=1
    for i in range(n):
        # evaluate at alpha^{-i} = alpha^{n - i}
        x = gf.exp[(n - i) % n]
        val = 0
        xp = 1
        for coeff in sigma:
            if coeff != 0:
                val ^= gf.mul(coeff, xp)
            xp = gf.mul(xp, x)
        if val == 0:
            # error at position i (LSB pos 0), map to index in codeword (MSB..LSB)
            idx = len(cw) - 1 - i
            if 0 <= idx < len(cw):
                cw[idx] ^= 1
                error_positions.append(idx)
    return cw, error_positions

# %% [markdown]
# بخش ۷ — دیکودر یکپارچه

# %%

def bch_decode(received: List[int], t: int, gf: GF2m) -> Tuple[List[int], bool, List[int]]:
    synd = compute_syndromes(received, t, gf)
    if all(s == 0 for s in synd):
        return received, True, []
    sigma = berlekamp_massey(synd, gf)
    corrected, errs = chien_search_and_correct(received, sigma, gf)
    # بررسی دوباره
    synd2 = compute_syndromes(corrected, t, gf)
    success = all(s == 0 for s in synd2)
    return corrected, success, errs

# %% [markdown]
# بخش ۸ — ابزار شبیه‌سازی کانال BSC و ارزیابی عملکرد

# %%

def flip_bits(codeword: List[int], p: float) -> Tuple[List[int], List[int]]:
    n = len(codeword)
    rx = codeword.copy()
    errs = []
    for i in range(n):
        if random.random() < p:
            rx[i] ^= 1
            errs.append(i)
    return rx, errs


def simulate_bch(m: int, t: int, prim_poly: int, trials: int = 1000, p_bit: float = 0.01, seed: int = None):
    if seed is not None:
        random.seed(seed)
        np.random.seed(seed)
    gf = GF2m(m, prim_poly)
    n = gf.n
    g = bch_generator_poly(m, t, prim_poly)
    k = n - (len(g) - 1)
    total_frame_err = 0
    total_bit_err = 0
    total_bits = trials * k
    for _ in range(trials):
        msg = [random.randint(0, 1) for _ in range(k)]
        cw = systematic_bch_encode(msg, g)
        rx, true_errs = flip_bits(cw, p_bit)
        decoded, success, found_errs = bch_decode(rx, t, gf)
        # compare decoded message
        dec_msg = decoded[:k]
        if dec_msg != msg:
            total_frame_err += 1
            # count bit errors in message portion
            for i in range(k):
                if dec_msg[i] != msg[i]:
                    total_bit_err += 1
    fer = total_frame_err / trials
    ber = total_bit_err / total_bits
    return {
        'm': m, 'n': n, 't': t, 'k': k, 'g': g,
        'trials': trials, 'p_bit': p_bit,
        'FER': fer, 'BER': ber
    }

# %% [markdown]
# بخش ۹ — مثال استفاده و نمایش نتایج

# %%
if __name__ == '__main__':
    # پارامترهای مثال — می‌توانید آن‌ها را تغییر دهید
    m = 4
    prim_poly = 0b10011  # x^4 + x + 1
    t = 2
    print(f"GF(2^{m}) with n = { (1<<m) - 1 }")
    g = bch_generator_poly(m, t, prim_poly)
    print("Generator polynomial g(x) (binary coeffs MSB..LSB):", g)
    k = ( (1<<m) - 1 ) - (len(g) - 1)
    print(f"Estimated k = {k} (information bits)")

    # یک قاب نمونه
    msg = [random.randint(0,1) for _ in range(k)]
    cw = systematic_bch_encode(msg, g)
    print("Message:", msg)
    print("Codeword:", cw)

    # وارد کردن خطا تا t بیت
    err_pos = random.sample(range(len(cw)), t)
    rx = cw.copy()
    for p in err_pos:
        rx[p] ^= 1
    print("Introduced errors at positions:", err_pos)

    decoded, success, found_errs = bch_decode(rx, t, GF2m(m, prim_poly))
    print("Decoding success:", success)
    print("Found error positions (indexes in codeword MSB..LSB):", found_errs)
    if success:
        print("Decoded message equals original?", decoded[:k] == msg)

    # شبیه‌سازی آماری کوتاه
    res = simulate_bch(m, t, prim_poly, trials=500, p_bit=0.02, seed=42)
    print("Simulation summary:")
    for k2,v in res.items():
        print(k2, ":", v)

# %% [markdown]
# توضیحات تکمیلی (پایین فایل)
#
# ۱) ساخت میدان GF(2^m):
#    - عناصر به صورت اعداد صحیح 0..2^m-1 نمایش داده می‌شوند؛ جدول exp/log برای ضرب سریع ساخته می‌شود.
#
# ۲) چندجمله‌ای مینیمال:
#    - برای هر ریشه alpha^i، ابتدا کلاس کونژوگاسیون (i, 2i mod n, 4i mod n, ...) محاسبه می‌شود.
#    - چندجمله‌ای مینیمال به‌صورت ضرب عوامل (x + alpha^{e}) در حوزهٔ ضرایب میدان محاسبه می‌شود.
#    - برای BCH دودویی، ضرایب نهایی چندجمله‌ای مینیمال باید در GF(2) باشند (0 یا 1).
#
# ۳) ساخت چندجمله‌ای مولد:
#    - مولد g(x) حاصل‌ضرب چندجمله‌ای‌های مینیمال برای ریشه‌های alpha^1 ... alpha^{2t} (فقط جذرهای جدید هر بار).
#
# ۴) انکُدینگ سیستماتیک:
#    - پیام k بیتی را به صورت پُر شده به درجهٔ deg(g) منتقل می‌کنیم، با g(x) تقسیم می‌کنیم و باقیمانده را انتها می‌چسبانیم.
#
# ۵) دیکُدینگ:
#    - سیندروم‌ها s_j = r(alpha^j) محاسبه می‌شوند (j=1..2t).
#    - الگوریتم Berlekamp–Massey بر روی این سیندروم‌ها اجرا می‌شود تا چندجمله‌ای مکان‌یاب خطا sigma(x) به‌دست آید.
#    - Chien search برای یافتن ریشه‌های sigma(x) (مکانهای خطا) مقداردهی می‌شود.
#    - برای کدهای باینری، مقدار خطا (error magnitude) برابر 1 است؛ لذا با یافتن مکان‌ها بیت‌ها را معکوس می‌کنیم.
#
# ۶) محدودیت‌ها و نکات عملی:
#    - این پیاده‌سازی آموزشی روی m‌های کوچکتر (m<=6 یا m<=7) مناسب است؛ برای m بزرگ‌تر لازم است
#      بهینه‌سازی‌های عددی و استفاده از ساختارهای دادهٔ سریع‌تر بکار رود.
#    - در برخی پیکربندی‌ها و یا اگر چندجمله‌ای مقدم انتخابی صحیح نباشد، ممکن است minimal_polynomial
#      تولید شود که ضرایب آن به‌صورت صحیح در GF(2) نباشند — برای پارامترهای متداول این مسئله رخ نمی‌دهد.
#    - برای استفاده تجاری یا پردازش‌های بزرگ توصیه می‌شود از کتابخانه‌های بهینه مثل `galois` یا
#      پیاده‌سازی‌های کامپایل شده استفاده کنید.
#
