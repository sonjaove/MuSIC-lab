# Project notes

## Dynamical Equations (19/11/25)

The event-triggered dynamics discussed on 19/11/25 were written in agent form as:

$$
\begin{cases}
\dot{x} = -\nabla f - Lx - Lz - L(\hat{x}-x) - L(\hat{z}-z),\\
\dot{z} = Lz + L(\hat{z}-z)
\end{cases}
$$

$$
\Rightarrow\quad
\begin{cases}
\dot{x} = -\nabla f - Lx - Lz - L e_x - L e_z,\\
\dot{z} = Lz + L e_z
\end{cases}
$$


<!-- In compact Laplacian form, these were expressed as:

$$
\dot{x} = -\nabla f(x) - L\hat{x} - L\hat{z},
$$

$$
\dot{z} = L\hat{z}.
$$

After expanding the sampled values  
$\hat{x} = x + e_x,\ \ \hat{z} = z + e_z$  
  
we obtained:

$$
\dot{x} = -\nabla f(x) - Lx - Lz - L e_x - L e_z,
$$

$$
\dot{z} = Lz + L e_z.
$$ -->

These are the equations exactly as written on the board on 19/11/25.



## Correction in Equations (10/12/25)

It was later identified that the expression for $\dot{z}$ was incorrect for the intended purpose.

The corrected form is [(has already been studied)](https://arxiv.org/abs/1407.7921) :

$$\quad
\begin{cases}
\dot{x} = -\nabla f - Lx - Lz - L e_x - L e_z,\\
\dot{z} = Lx + L e_x
\end{cases}
$$

Here:

- $z$ now behaves as an **integrator of disagreement in $x$**,  
- the term $-z$ in $\dot{x}$ provides the required **integral correction**,  
- the system now structurally matches the role of $v$ in the original algorithm.

This correction restores the intended convergence behaviour.

## Results
1. ![pairwise dis](pairwise_disagreement.png)
2. ![trajectroies](<trejactories and component-wise convergence.png>)
- the only caveat is that the there is some oscillation when agents converge

## recent advancments in this field :
> google search: "recent advancements in the field of Distributed convex optimization via continuous-time coordination algorithms with discrete-time communication"
1. using PI controller instead of event-triggered contorl.
2. subgradient methods
3. [Distributed Predefined-time Zero-gradient-sum Optimization for Multi-agent Systems: From Continuous-time to Event-triggered Communication (recent most 2025)](https://link.springer.com/article/10.1007/s12555-023-0140-1)
