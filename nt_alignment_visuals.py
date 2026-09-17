#!/usr/bin/env python3
import sys
from Bio import AlignIO
import matplotlib.pyplot as plt
import numpy as np

# --- Nature Formatting Configuration ---
INPUT_FASTA = "path/to/alignment.fasta"  # Replace with your input FASTA file path
OUTPUT_PDF = "visual.pdf"
START_POS = 1  # 1-based start
END_POS = 1000  # 1-based end - END_POS should be the last position you want to include in the visualisation
MAX_BLOCKS = 4  # Exactly 4 stacked blocks - change this depending how long your gene is e.g. for a 1.5kb gene, 4 blocks of 400bp each is ideal for Nature figures

# Typography rules (Min 5pt, Max 7pt)
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = [
    "Helvetica",
    "DejaVu Sans",
    "Liberation Sans",
    "sans-serif",
]
plt.rcParams["pdf.fonttype"] = 42
FONT_SIZE_LABEL = 6.5  # pt (Sequence names)
FONT_SIZE_NUM = 5.5  # pt (Coordinates & ticks)

# Color Palette (Colorblind-friendly, IBM palette)
PALETTE = {
    "A": [0.392, 0.561, 1.000],  # Blue    (#648FFF)
    "T": [1.000, 0.690, 0.000],  # Gold/Yellow (#FFB000)
    "C": [0.863, 0.149, 0.498],  # Magenta (#DC267F)
    "G": [0.996, 0.380, 0.000],  # Orange  (#FE6100)
    "-": [1.000, 1.000, 1.000],  # White (Gaps)
    "N": [0.600, 0.600, 0.600],  # Grey
}
DEFAULT_COLOR = [1.000, 1.000, 1.000]

def render_alignment():
  alignment = AlignIO.read(INPUT_FASTA, "fasta")
  num_seqs = len(alignment)
  aln_len = alignment.get_alignment_length()

  start = max(0, START_POS - 1)
  end = min(aln_len, END_POS)
  sub_len = end - start

  # Distribute across a maximum of 4 blocks
  res_per_block = int(np.ceil(sub_len / MAX_BLOCKS))
  num_blocks = int(np.ceil(sub_len / res_per_block))

  # Dimensions in inches for Nature double-column (183 mm = 7.205 inches)
  fig_width_in = 183 / 25.4

  # Height estimation: tightly scaled without legend margin
  block_height_in = (num_seqs * 0.08 + 0.30) * (72 / 25.4) / 72
  fig_height_in = max(2.2, min(7.0, block_height_in * num_blocks + 0.3))

  fig, axes = plt.subplots(
      num_blocks, 1, figsize=(fig_width_in, fig_height_in), squeeze=False
  )

  for block_idx in range(num_blocks):
    ax = axes[block_idx, 0]
    b_start = start + block_idx * res_per_block
    b_end = min(start + (block_idx + 1) * res_per_block, end)
    cur_len = b_end - b_start

    # Build RGB matrix for the block
    grid_rgb = np.zeros((num_seqs, cur_len, 3))
    for seq_idx, record in enumerate(alignment):
      sub_seq = str(record.seq[b_start:b_end]).upper()
      for x, char in enumerate(sub_seq):
        grid_rgb[seq_idx, x] = PALETTE.get(char, DEFAULT_COLOR)

    # Render colored heatmap matrix without letters
    ax.imshow(
        grid_rgb,
        aspect="auto",
        interpolation="nearest",
        extent=[-0.5, cur_len - 0.5, num_seqs - 0.5, -0.5],
    )

    # Sequence ID labels on Y-axis
    ax.set_yticks(range(num_seqs))
    ax.set_yticklabels(
        [r.id for r in alignment], fontsize=FONT_SIZE_LABEL, fontweight="normal"
    )
    ax.tick_params(
        axis="y",
        length=0,
        pad=3,
        labelsize=FONT_SIZE_LABEL,
        labelcolor="#222222",
    )

    # Top position ruler ticks
    tick_step = 50 if cur_len > 150 else 20
    x_ticks = np.arange(0, cur_len, tick_step)
    x_labels = [str(b_start + x + 1) for x in x_ticks]

    ax.set_xticks(x_ticks)
    ax.set_xticklabels(x_labels, fontsize=FONT_SIZE_NUM, color="#444444")
    ax.xaxis.tick_top()
    ax.tick_params(axis="x", length=2, pad=1, direction="out", color="#888888")

    # Clean styling
    for spine in ax.spines.values():
      spine.set_edgecolor("#CCCCCC")
      spine.set_linewidth(0.5)

    ax.set_xlim(-0.5, cur_len - 0.5)

  # Adjusted margins to reclaim bottom space previously occupied by the key
  plt.subplots_adjust(
      left=0.18, right=0.96, top=0.94, bottom=0.04, hspace=0.50
  )
  plt.savefig(OUTPUT_PDF, format="pdf", bbox_inches="tight")
  print(f"Publication-ready Nature PDF saved to: {OUTPUT_PDF}")


if __name__ == "__main__":
  render_alignment()
  