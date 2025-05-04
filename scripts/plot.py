import pandas as pd
import matplotlib.pyplot as plt
import os

DATA_DIR = "data"
FIGS_DIR = "figs"

if not os.path.exists(FIGS_DIR):
    os.makedirs(FIGS_DIR)

SUMMARY_FILE = os.path.join(DATA_DIR, "summary.csv")
RESIDUALS_FILE = os.path.join(DATA_DIR, "residuals.csv")
COMPARISON_FILE = os.path.join(DATA_DIR, "comparison.csv")
TARGET_RESIDUAL = 1e-7

df_summary = pd.read_csv(SUMMARY_FILE)
df_residuals = pd.read_csv(RESIDUALS_FILE)

fig1, (ax1a, ax1b) = plt.subplots(1, 2, figsize=(12, 5))
fig1.suptitle(r"Performance vs. Number of Levels ($N=128$)")

ax1a.plot(df_summary["Levels"], df_summary["Cycles"], marker="o", linestyle="-")
ax1a.set_xlabel("lmax")
ax1a.set_ylabel("Cycles")
ax1a.set_title("Cycles vs. Levels")
ax1a.grid(True, which="both", linestyle="--", linewidth=0.5)
if not df_summary.empty:
    ax1a.set_xticks(df_summary["Levels"])

ax1b.plot(
    df_summary["Levels"],
    df_summary["Runtime"],
    marker="o",
    linestyle="-",
    color="r",
)
ax1b.set_xlabel("lmax")
ax1b.set_ylabel("Runtime (s)")
ax1b.set_title("Runtime vs. Levels")
ax1b.grid(True, which="both", linestyle="--", linewidth=0.5)
if not df_summary.empty:
    ax1b.set_xticks(df_summary["Levels"])

fig1.tight_layout(rect=[0, 0.03, 1, 0.95])
plt.savefig(os.path.join(FIGS_DIR, "q1_perf-vs-lmax.png"))
plt.close(fig1)

available_lmax_values = sorted(df_residuals["Levels"].unique())

for lmax_val in available_lmax_values:
    df_lmax_conv = df_residuals[df_residuals["Levels"] == lmax_val].copy()
    df_lmax_conv.sort_values(by="CycleIteration", inplace=True)

    if not df_lmax_conv.empty:
        plt.figure(figsize=(8, 6))

        plt.semilogy(
            df_lmax_conv["CycleIteration"],
            df_lmax_conv["ResidualNorm"],
            marker=".",
            markersize=3,
            linestyle="-",
            label=f"lmax = {lmax_val}",
        )

        plt.axhline(
            y=TARGET_RESIDUAL,
            color="grey",
            linestyle=":",
            linewidth=1.5,
            label="Target Residual",
        )
        plt.xlabel("Cycles")
        plt.ylabel("2-norm")
        plt.title(rf"Convergence History for $N=128$, lmax = {lmax_val}")

        max_cycles = df_lmax_conv["CycleIteration"].max()
        min_resid = df_lmax_conv["ResidualNorm"].min()
        initial_resid_series = df_lmax_conv[df_lmax_conv["CycleIteration"] == 0][
            "ResidualNorm"
        ]
        max_resid = (
            initial_resid_series.iloc[0]
            if not initial_resid_series.empty
            else df_lmax_conv["ResidualNorm"].max()
        )

        plt.xlim(left=-0.02 * max_cycles, right=max_cycles * 1.05)
        plt.ylim(
            bottom=max(TARGET_RESIDUAL / 10, min_resid / 10),
            top=max_resid * 1.5,
        )

        plt.legend()
        plt.grid(True, which="both", linestyle="--", linewidth=0.5)
        plt.tight_layout()

        plot_filename = f"q1_residuals-lmax{lmax_val}.png"
        plt.savefig(os.path.join(FIGS_DIR, plot_filename))
        plt.close()

df_comparison = pd.read_csv(COMPARISON_FILE)

df_comparison = df_comparison[df_comparison["Cycles"] > 0].copy()

df_2level = df_comparison[df_comparison["Config"] == "2-level"].sort_values("N")
df_maxlevel = df_comparison[df_comparison["Config"] == "max-level"].sort_values("N")

fig2, (ax2a, ax2b) = plt.subplots(1, 2, figsize=(14, 6))
fig2.suptitle("2-level vs. max-level")

ax2a.plot(
    df_2level["N"],
    df_2level["Cycles"],
    marker="o",
    linestyle="-",
    label="2-Level",
)
ax2a.plot(
    df_maxlevel["N"],
    df_maxlevel["Cycles"],
    marker="s",
    linestyle="--",
    label="Max-Level (N=8 Coarsest)",
)
ax2a.set_xlabel("Grid Size (N_interior)")
ax2a.set_ylabel("Cycles")
ax2a.set_title("Cycles vs. Grid Size")
ax2a.set_xscale("log", base=2)
ax2a.set_xticks(df_2level["N"])
ax2a.set_xticklabels(df_2level["N"])
ax2a.legend()
ax2a.grid(True, which="both", linestyle="--", linewidth=0.5)

ax2b.plot(
    df_2level["N"],
    df_2level["Runtime"],
    marker="o",
    linestyle="-",
    label="2-Level",
)
ax2b.plot(
    df_maxlevel["N"],
    df_maxlevel["Runtime"],
    marker="s",
    linestyle="--",
    label="Max-Level (N=8 Coarsest)",
)
ax2b.set_xlabel("Grid Size (N_interior)")
ax2b.set_ylabel("Runtime (s)")
ax2b.set_title("Runtime vs. Grid Size")
ax2b.set_xscale("log", base=2)
ax2b.set_yscale("log")
ax2b.set_xticks(df_2level["N"])
ax2b.set_xticklabels(df_2level["N"])
ax2b.legend()
ax2b.grid(True, which="both", linestyle="--", linewidth=0.5)

fig2.tight_layout(rect=[0, 0.03, 1, 0.95])
plt.savefig(os.path.join(FIGS_DIR, "q2_comparison.png"))
plt.close(fig2)
