---
tags:
  - computational-model
---
# SAC: Spacing effects investigations

The Source of Activation Confusion (SAC) is a mechanistic model of declarative memory. It was developed originaly by Reder and colleagues to explain cognitive illusions (Reder & Gordon, 1997) and feeling-of-knowing (Reder & Schunn, 1996). Since then it has accounted for many key findings in human memory, such as effects of mirror word frequency (Reder et al., 2000), list length (Cary & Reder, 2003), and list strength (Diana & Reder, 2005), as well as how prior knowledge affects learning (Reder et al, 2007; Popov & Reder, 2020). Most recently, SAC has helped us understand how people use their limited cognitive resources to prioritize relevant information (Popov et al, 2019) and to maximize future rewards (Ma, Popov, & Zhang, 2026). 

Despite these successes, SAC has so far been unable to explain one of the oldest and most robust effects in human memory - the spacing effect. SAC shares roots with ACT-R’s declarative module, and while the two models differ in key respects, until recently SAC built upon ACT-R’s core learning equations. These equations do not produce spacing effects in their standard form (Pavlik & Anderson, 2005; Anderson et al, 2023). 

However, in 2020, building on earlier insights about the role of familiarity on learning (i.e. Diana & Reder, 2006; Reder et al, 2007), we introduced many changes to SAC’s dynamics, including what seemed like a minor adjustment to the equation governing memory strength. Unexpectedly, the revised SAC model reproduces many key spacing results.

Mechanistic models like SAC have many moving parts and it is not always clear which components and processes are responsible for a model’s behavior. Here we strip down the 2020 SAC model to the essential equation responsible for its novel spacing behavior. We examine its theoretical relations to existing models of spacing effects, and investigate its properties.

## Core model

SAC implements a spreading activation theory in which semantic, episodic and contextual memory traces are represented as localist nodes in a network. Each node (and the links between them) has a continuous scalar memory strength value. As in ACT-R, this base-level strength reflects the history of use - it increases through practice and decays with time. 

Specifically, each study event creates a strength increment whose size depends on the current total strength, and every increment decays independently. 

Let $0 =t_1 < t_2 < \cdots < t_n < \cdots$ be the absolute times at which a repeated event occurs. With $B(t)$ as the base-level strength at $t$, $u_n$ as the increment created by the $n$\-th presentation, and a learning rate $\delta \in (0,1]$, we have

$$
u_n=\delta(1-B(t_n)), \quad B(t)=\sum_{k: \,t_k<t} u_k f(t-t_k), \quad B(0) = 0
$$

Here, $f$ is the increment decay function, which is a shifted power law:

$$
f(t) = \left(1 + t \right)^{-d}
$$

As shown in Derive SAC’s base-level equation full form, we get the following direct form for the base-level strength at a time of testing $t > t_n$ after $n$ repetitions:

$$
B_n(t) = \sum_{j=1}^{n}(-1)^{j+1}\delta^j \sum_{1 \le k_1 < k_2 < \cdots < k_j \le n}  f(t-t_{k_j})  f(t_{k_j}-t_{k_{j-1}}) \cdots f(t_{k_2}-t_{k_1})
$$

Examples:

$$
\begin{aligned}
B_1(t) &= \delta f(t -t_1) \\
B_2(t) &= \delta(f(t-t_2) +f(t-t_1)) \\ &- \delta^2 f(t-t_2)f(t_2-t_1) \\
B_3(t) &= \delta \left[f(t-t_3) +f(t-t_2) + f(t-t_1) \right] \\ &- \delta^2 \left[f(t-t_3) f(t_3-t_2) + f(t-t_3)f(t_3 - t_1) + f(t-t_2)(t_2-t_1)\right] \\ &+ \delta^3f(t-t_3)f(t_3-t_2)f(t_2 - t_1)
\end{aligned}
$$

## Optimal spacing for two study events

Assume a standard setup where we have two study events followed by a test. Let $x$ be the time lag between the two study events and $y$ be the retention interval. Classic results show that performance is non-monotonic - the optimal study lag increases as the retention interval increases.

In SAC, the strength is:

$$
B(x, y) = \delta f(y) +\delta f(x+y) - \delta^2 f(x)f(y)
$$

This equation has an interior optimum $x^*$ if it satisfies:

$$
f'(x^*+y) = \delta f'(x^*) f(y), \quad f''(x^*+y) < \delta f''(x^*) f(y)
$$

For the shifted power rule we get for the first condition:

$$
(1+x^*+y)^{-d-1}= \delta(1+x^*)^{-d-1}(1+y)^{-d}
$$

$$
1+x^*+y =\delta^{-1/(d+1)}(1+x^*)(1+y)^{d/(d+1)}
$$

$$
y = (1+x^*)(\delta^{-1/(d+1)}(1+y)^{d/(d+1)}-1)
$$

$$
x^* = \frac{y}{\delta^{-1/(d+1)}(1+y)^{d/(d+1)}-1} - 1
$$

This must be positive, which holds for

$$
y > \frac{1}{\delta}-1
$$

For large retention intervals $y \gg 1$

$$
x^*
\sim
\left(\delta y\right)^{1/(d+1)}.
$$

hence

$$
\log x^* \sim \frac{1}{1+d} \log y + \mathrm{constant}
$$





