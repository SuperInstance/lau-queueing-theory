# lau-queueing-theory

A pure-Rust library for the mathematical analysis of waiting-line (queueing) systems — from the classic M/M/1 and M/M/c models through Jackson networks, Erlang loss/delay formulas, birth-death processes, and priority-queue scheduling.

**72 tests** · `nalgebra` + `serde` · MIT license

---

## What This Does

`lau-queueing-theory` provides closed-form and algorithmic solutions for the performance metrics of queueing systems: expected queue length, waiting time, server utilisation, blocking probability, throughput, and more.

You define a system by its **Kendall notation** (A/S/c/K/N/D), call the appropriate module, and get back typed structs with every metric you'd compute by hand in an operations-research course. No simulation, no Monte Carlo — exact formulas where they exist, numerically stable recursions where they don't.

---

## Key Idea

Every queueing system reduces to a few parameters: arrival rate λ, service rate μ, number of servers *c*, system capacity *K*, population size *N*, and discipline (FIFO, LIFO, priority). This crate encodes the **standard analytical results** from Gross & Harris, Kleinrock, and Bolch et al. as type-safe Rust functions:

| Module | System | What you get |
|---|---|---|
| `kendall` | Kendall notation parser/builder | Structured `KendallNotation` |
| `mm1` | M/M/1 | L, L_q, W, W_q, ρ, P_n, P(N>n) |
| `mmc` | M/M/c (multi-server) | Same metrics, Erlang-C formula |
| `mg1` | M/G/1 (Pollaczek–Khinchine) | L, L_q, W, W_q from variance of service time |
| `erlang` | Erlang-B (loss), Erlang-C (delay) | Blocking probability, waiting probability |
| `birth_death` | General birth-death | Steady-state probabilities via recurrence |
| `networks` | Jackson networks | Throughput, queue lengths, visit ratios |
| `priority` | M/M/1 & M/M/c with priorities | Per-class metrics (preemptive & non-preemptive) |
| `agent_scheduling` | Agent/server scheduling | Optimal staffing, cost models |

---

## Install

```toml
[dependencies]
lau-queueing-theory = "0.1"
```

Or:

```sh
cargo add lau-queueing-theory
```

---

## Quick Start

### M/M/1 Queue

```rust
use lau_queueing_theory::mm1::MM1Result;

let result = MM1Result::from_lambda_mu(5.0, 8.0);
// λ = 5 arrivals/sec, μ = 8 services/sec

println!("Utilisation ρ = {:.3}", result.rho);           // 0.625
println!("Mean queue length L_q = {:.3}", result.lq);    // ~1.042
println!("Mean waiting time W_q = {:.3} s", result.wq);  // ~0.208 s
println!("P(empty) = {:.3}", result.pn(0));               // 0.375
```

### Erlang-B: How many circuits?

```rust
use lau_queueing_theory::erlang::erlang_b;

let traffic = 12.0;  // 12 Erlangs of offered load
let servers = 20;    // try 20 circuits
let blocking = erlang_b(traffic, servers);
println!("Blocking probability with {} servers: {:.4}", servers, blocking);
```

### Jackson Network

```rust
use lau_queueing_theory::networks::{JacksonNetwork, Node};

let mut net = JacksonNetwork::new();
net.add_node(Node::new(0, 1.0, 2.0));  // node 0: λ_ext=1, μ=2
net.add_node(Node::new(1, 0.0, 3.0));  // node 1: no ext. arrivals, μ=3
net.set_routing(0, 1, 0.6);            // 60% from node 0 → node 1
net.set_routing(1, 0, 0.3);            // 30% from node 1 → node 0

let solution = net.solve();
for (i, node_result) in solution.iter().enumerate() {
    println!("Node {}: λ_eff={:.3}, L={:.3}, W={:.3}",
        i, node_result.lambda_eff, node_result.l, node_result.w);
}
```

---

## API Reference

### `kendall` — Kendall Notation

```rust
pub struct KendallNotation {
    pub arrival: ArrivalProcess,   // M, D, Ek, G, ...
    pub service: ServiceTime,      // M, D, Ek, G, ...
    pub servers: usize,            // c
    pub capacity: Option<usize>,   // K (default ∞)
    pub population: Option<usize>, // N (default ∞)
    pub discipline: Discipline,    // FIFO, LIFO, SIRO, PRI
}
```

Parse from string: `KendallNotation::parse("M/M/c/K/N/FIFO")`.

### `mm1` — M/M/1 Queue

| Method / Field | Description |
|---|---|
| `MM1Result::from_lambda_mu(λ, μ)` | Construct from arrival & service rates |
| `.rho` | Server utilisation ρ = λ/μ |
| `.l` | Mean number in system L = λ/(μ−λ) |
| `.lq` | Mean queue length L_q = λ²/(μ(μ−λ)) |
| `.w` | Mean time in system W = 1/(μ−λ) |
| `.wq` | Mean waiting time W_q = λ/(μ(μ−λ)) |
| `.pn(n)` | P(N=n) = (1−ρ)ρⁿ |
| `.p_n_gt(k)` | P(N>k) = ρ^{k+1} |

### `mmc` — M/M/c Queue

| Method / Field | Description |
|---|---|
| `MMCResult::from_lambda_mu_c(λ, μ, c)` | Construct |
| `.rho`, `.l`, `.lq`, `.w`, `.wq` | Standard metrics |
| `.erlang_c()` | P(waiting) = Erlang-C formula |
| `.pn(n)` | Steady-state P(N=n) |

### `mg1` — M/G/1 Queue (Pollaczek–Khinchine)

```rust
pub struct MG1Result {
    pub rho: f64, pub l: f64, pub lq: f64,
    pub w: f64, pub wq: f64,
}
```

Construct via `MG1Result::from_params(λ, E[S], Var[S])` where E[S] is mean service time and Var[S] is its variance.

**Pollaczek–Khinchine formula:** L_q = λ²(Var[S] + E[S]²) / (2(1 − ρ))

### `erlang` — Erlang Formulas

| Function | Signature | Description |
|---|---|---|
| `erlang_b(a, c)` | `(f64, usize) → f64` | Erlang-B: blocking probability for M/M/c/c (loss system) |
| `erlang_c(a, c)` | `(f64, usize) → f64` | Erlang-C: delay probability for M/M/c |

### `birth_death` — General Birth-Death Processes

```rust
pub struct BirthDeath {
    pub lambda: Vec<f64>,  // birth rates λ_0, λ_1, ..., λ_{n-1}
    pub mu: Vec<f64>,      // death rates μ_1, μ_2, ..., μ_n
}

impl BirthDeath {
    pub fn steady_state(&self) -> Vec<f64>;
    pub fn mean_population(&self) -> f64;
    pub fn mean_rate(&self) -> f64;
}
```

Steady-state computed via the standard recurrence: π_n = (λ₀λ₁…λ_{n-1})/(μ₁μ₂…μ_n) · π_0, with π_0 determined by normalisation.

### `networks` — Jackson Networks

| Method | Description |
|---|---|
| `JacksonNetwork::new()` | Create empty network |
| `.add_node(Node)` | Add a service station |
| `.set_routing(from, to, p)` | Set routing probability |
| `.solve()` | Solve traffic equations → per-node metrics |

Solves the traffic flow equations **λ_i = λ_ext,i + Σ_j λ_j · p_{ji}** via matrix inversion, then treats each node as an independent M/M/c queue.

### `priority` — Priority Queues

```rust
pub struct PriorityQueue {
    pub classes: Vec<PriorityClass>,
}

pub struct PriorityClass {
    pub lambda: f64, pub mu: f64,
    pub lq: f64, pub wq: f64, pub w: f64,
}
```

Supports both **preemptive** and **non-preemptive** resume disciplines. Higher-priority classes have lower index.

### `agent_scheduling` — Agent/Server Scheduling

```rust
pub struct AgentSchedule {
    pub num_agents: usize,
    pub cost_per_agent: f64,
    pub waiting_cost_rate: f64,
    pub total_cost: f64,
}

pub fn optimal_staffing_mms(lambda: f64, mu: f64, s_max: usize,
    cost_agent: f64, cost_waiting: f64) -> AgentSchedule;
```

Finds the number of agents *s* that minimises total cost = s·cost_agent + L_q·cost_waiting.

---

## How It Works

### Steady-State Probability Computation

For M/M/c systems, steady-state probabilities are computed using the recursive formula:

- π_0 is found via normalisation
- π_n for n ≤ c uses the standard term (λ/μ)ⁿ/n!
- π_n for n > c uses (λ/μ)ⁿ/(c!·c^{n−c})

Numerical stability is maintained by computing ratios rather than large factorials.

### Erlang-B Recursion

The Erlang-B formula is computed via the efficient recursion:

**B(0, a) = 0**, **B(c, a) = a·B(c−1, a) / (1 + a·B(c−1, a))**

This avoids overflow from computing large factorials directly.

### Jackson Network Traffic Equations

For a network with routing matrix P and external arrivals **λ_ext**, the effective arrival rates satisfy:

**Λ = λ_ext + Λ · P**

Solved as **Λ = λ_ext · (I − P)⁻¹** using `nalgebra` matrix inversion.

---

## The Math

### Little's Law

For any queueing system in steady state: **L = λW**

Mean number in system = arrival rate × mean time in system. This universal result underpins every metric in the crate.

### Pollaczek–Khinchine (M/G/1)

For an M/G/1 queue with arrival rate λ, service time mean E[S] and variance Var[S]:

- ρ = λE[S]
- L_q = (λ² · Var[S] + ρ²) / (2(1 − ρ))
- W_q = L_q / λ

### Erlang-B (M/M/c/c Loss System)

Probability all *c* servers are busy (call is lost):

**B(c, a) = (a^c / c!) / Σ_{k=0}^{c} a^k / k!**

where a = λ/μ is the offered traffic in Erlangs.

### Jackson's Theorem

For an open Jackson network of M/M/c nodes with probabilistic routing, each node behaves as an independent M/M/c queue with effective arrival rate given by the traffic equations. The network steady-state distribution is the product of individual node distributions.

---

## License

MIT
