# Project notes

## Dynamical Equations (19/11/25)

The event-triggered dynamics discussed on 19/11/25 were written in agent form as:

$$
\begin{cases}
\dot{x} = -\nabla f - Lx - Lz - L(\hat{x}-x) - L(\hat{z}-z),\\[6pt]
\dot{z} = Lz + L(\hat{z}-z)
\end{cases}
$$

$$
\Rightarrow\quad
\begin{cases}
\dot{x} = -\nabla f - Lx - Lz - L e_x - L e_z,\\[6pt]
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



## Correction in Equations

It was later identified that the expression for $\dot{z}$ was incorrect for the intended purpose.

The corrected form is:

$$
\begin{cases}
\dot{x} = -\nabla f - Lx - Lz - L e_x - L e_z, \\[6pt]
\dot{z} = Lx + L e_x
\end{cases}
$$
Here:

- $z$ now behaves as an **integrator of disagreement in $x$**,  
- the term $-z$ in $\dot{x}$ provides the required **integral correction**,  
- the system now structurally matches the role of $v$ in the original algorithm.

This correction restores the intended convergence behaviour.
- $z$ now behaves as an **integrator of disagreement in $x$**,  
- the term $-z$ in $\dot{x}$ provides the required **integral correction**,  
- the system now structurally matches the role of $v$ in the original algorithm.

This correction restores the intended convergence behaviour.

## Results
1. ![pairwise dis](pairwise_disagreement.png)
2. ![trajectroies](<trejactories and component-wise convergence.png>)
- the only caveat is that the there is some oscillation when agents converge

