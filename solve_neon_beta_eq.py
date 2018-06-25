from sympy import (symbols, exp, Ynm, Abs, I, pi, cos, sin, sqrt,
                   cancel, expand_func, legendre, solve, re, simplify, expand,
                   init_printing, pprint, latex)

init_printing()


# %% expand left term
phi, the, vphi = symbols('phi theta varphi', real=True)
eta_p, eta_f, eta_s, eta_d = symbols(
    'eta_p eta_f eta_s eta_d',
    real=True,
)
c_psp, c_pdp, c_fdp, c_sp, c_dp = symbols(
    'c_psp c_pdp c_fdp c_sp c_dp',
    positive=True,
)

c0sp = -sqrt(3)/3*c_sp
c0dp, c1dp = sqrt(6)/3*c_dp, sqrt(2)/2*c_dp
c0psp = -sqrt(3)/3*c_psp
c0pdp, c1pdp = -2*sqrt(15)/15*c_pdp, -sqrt(15)/10*c_pdp
c0fdp, c1fdp = sqrt(10)/5*c_fdp, 2*sqrt(15)/15*c_fdp
waves_m = {
    -1: (c1pdp * exp(eta_p * I) * Ynm(1, -1, the, vphi) +
         c1fdp * exp(eta_f * I) * Ynm(3, -1, the, vphi) +
         c1dp * exp(eta_d * I + phi * I) * Ynm(2, -1, the, vphi)),
    0: (c0psp * exp(eta_p * I) * Ynm(1, 0, the, vphi) +
        c0pdp * exp(eta_p * I) * Ynm(1, 0, the, vphi) +
        c0fdp * exp(eta_f * I) * Ynm(3, 0, the, vphi) +
        c0sp * exp(eta_s * I + phi * I) * Ynm(0, 0, the, vphi) +
        c0dp * exp(eta_d * I + phi * I) * Ynm(2, 0, the, vphi)),
    1: (c1pdp * exp(eta_p * I) * Ynm(1, 1, the, vphi) +
        c1fdp * exp(eta_f * I) * Ynm(3, 1, the, vphi) +
        c1dp * exp(eta_d * I + phi * I) * Ynm(2, 1, the, vphi)),
}
yield_lft = (
    expand_func(Abs(waves_m[-1]) ** 2 + Abs(waves_m[0]) ** 2 + Abs(waves_m[1]) ** 2)
    .subs(sin(the) ** 2, 1 - cos(the) ** 2)
)
#選択するmを選択する
__expr = cancel(yield_lft)
term0_lft = __expr.subs(the, pi / 2)
__expr = cancel((__expr - term0_lft) / cos(the))
term1_lft = __expr.subs(the, pi / 2)
__expr = cancel((__expr - term1_lft) / cos(the))
term2_lft = __expr.subs(the, pi / 2)
__expr = cancel((__expr - term2_lft) / cos(the))
term3_lft = __expr.subs(the, pi / 2)
__expr = cancel((__expr - term3_lft) / cos(the))
term4_lft = __expr.subs(the, pi / 2)
__expr = cancel((__expr - term4_lft) / cos(the))
term5_lft = __expr.subs(the, pi / 2)
__expr = cancel((__expr - term5_lft) / cos(the))
term6_lft = __expr.subs(the, pi / 2)
__expr = cancel((__expr - term6_lft) / cos(the))
if __expr != 0:
    AssertionError('Valuable __expr is not zero!')


# %% expand right term
b0, b1, b2, b3, b4, b5, b6 = symbols(
    'b_0 b_1 b_2 b_3 b_4 b_5 b_6',
    real=True
)
yield_rgt = (b0 + b1 * legendre(1, cos(the)) + b2 * legendre(2, cos(the)) +
             b3 * legendre(3, cos(the)) + b4 * legendre(4, cos(the)) +
             b5 * legendre(5, cos(the)) + b6 * legendre(6, cos(the)))

__expr = cancel(yield_rgt)
term0_rgt = __expr.subs(the, pi / 2)
__expr = cancel((__expr - term0_rgt) / cos(the))
term1_rgt = __expr.subs(the, pi / 2)
__expr = cancel((__expr - term1_rgt) / cos(the))
term2_rgt = __expr.subs(the, pi / 2)
__expr = cancel((__expr - term2_rgt) / cos(the))
term3_rgt = __expr.subs(the, pi / 2)
__expr = cancel((__expr - term3_rgt) / cos(the))
term4_rgt = __expr.subs(the, pi / 2)
__expr = cancel((__expr - term4_rgt) / cos(the))
term5_rgt = __expr.subs(the, pi / 2)
__expr = cancel((__expr - term5_rgt) / cos(the))
term6_rgt = __expr.subs(the, pi / 2)
__expr = cancel((__expr - term6_rgt) / cos(the))
if __expr != 0:
    AssertionError('Valuable __expr is not zero!')


# %%
b6_cmpx = simplify(cancel(solve(term6_lft - term6_rgt, b6)[0]))
b6_real = simplify(re(expand(b6_cmpx)))
print('b6 =')
pprint(b6_real, use_unicode=False)
print(latex(b6_real))
print('')

b5_cmpx = simplify(cancel(solve(term5_lft - term5_rgt, b5)[0]))
b5_real = simplify(re(expand(b5_cmpx)))
print('b5 =')
pprint(b5_real, use_unicode=False)
print(latex(b5_real))
print('')

b4_cmpx = simplify(cancel(solve((term4_lft - term4_rgt)
                                .subs(b6, b6_cmpx), b4)[0]))
b4_real = simplify(re(expand(b4_cmpx)))
print('b4 =')
pprint(b4_real, use_unicode=False)
print(latex(b4_real))
print('')

b3_cmpx = simplify(cancel(solve((term3_lft - term3_rgt)
                                .subs(b5, b5_cmpx), b3)[0]))
b3_real = simplify(re(expand(b3_cmpx)))
print('b3 =')
pprint(b3_real, use_unicode=False)
print(latex(b3_real))
print('')

b2_cmpx = simplify(cancel(solve((term2_lft - term2_rgt)
                                .subs(b6, b6_cmpx)
                                .subs(b4, b4_cmpx), b2)[0]))
b2_real = simplify(re(expand(b2_cmpx)))
print('b2 =')
pprint(b2_real, use_unicode=False)
print(latex(b2_real))
print('')

b1_cmpx = simplify(cancel(solve((term1_lft - term1_rgt)
                                .subs(b5, b5_cmpx)
                                .subs(b3, b3_cmpx), b1)[0]))
b1_real = simplify(re(expand(b1_cmpx)))
print('b1 =')
pprint(b1_real, use_unicode=False)
print(latex(b1_real))
print('')

b0_cmpx = simplify(cancel(solve((term0_lft - term0_rgt)
                                .subs(b6, b6_cmpx)
                                .subs(b4, b4_cmpx)
                                .subs(b2, b2_cmpx), b0)[0]))
b0_real = simplify(re(expand(b0_cmpx)))
print('b0 =')
pprint(b0_real, use_unicode=False)
print(latex(b0_real))
print('')
