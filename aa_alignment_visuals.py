#!/usr/bin/env python3
import sys
import numpy as np
import matplotlib.pyplot as plt
from Bio import AlignIO

# --- Nature Formatting Configuration ---
INPUT_FASTA = "path/to/alignment.fasta"  # Replace with your input FASTA file path
OUTPUT_PDF  = "visual.pdf"
START_POS   = 1       # 1-based start residue
END_POS     = 483     # 1-based end residue - END_POS should be the last position you want to include in the visualisation
MAX_BLOCKS  = 4       # Max 4 stacked blocks - change this depending how long your gene is, e.g. a 100aa gene could be visualised in 1 block.

# Nature typography rules (Min 5pt, Max 7pt)
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Helvetica', 'DejaVu Sans', 'Liberation Sans', 'sans-serif']
plt.rcParams['pdf.fonttype'] = 42
FONT_SIZE_LABEL = 6.5
FONT_SIZE_NUM   = 5.5

# Unique distinct color for every individual amino acid
AA_PALETTE = {
    'A': [0.55, 0.70, 0.00],  # Olive green
    'C': [1.00, 0.85, 0.00],  # Bright Gold / Yellow
    'D': [0.90, 0.10, 0.10],  # Vivid Red
    'E': [0.75, 0.00, 0.20],  # Crimson / Dark Red
    'F': [0.20, 0.40, 0.80],  # Royal Blue
    'G': [0.90, 0.50, 0.10],  # Orange
    'H': [0.50, 0.30, 0.70],  # Slate Purple
    'I': [0.10, 0.70, 0.80],  # Cyan / Cerulean
    'K': [0.60, 0.10, 0.80],  # Bright Purple
    'L': [0.15, 0.50, 0.75],  # Medium Blue
    'M': [0.85, 0.70, 0.10],  # Ochre / Mustard
    'N': [0.40, 0.80, 0.60],  # Mint Green
    'P': [0.60, 0.60, 0.60],  # Medium Gray
    'Q': [0.10, 0.65, 0.50],  # Sea green
    'R': [0.35, 0.10, 0.65],  # Deep Violet
    'S': [0.95, 0.65, 0.75],  # Rose / Pink
    'T': [0.85, 0.45, 0.55],  # Dusty Coral
    'V': [0.30, 0.75, 0.40],  # Leaf Green
    'W': [0.20, 0.20, 0.60],  # Dark Navy
    'Y': [0.45, 0.60, 0.85],  # Steel Blue
    '-': [1.00, 1.00, 1.00],  # White (Gaps)
    'X': [0.85, 0.85, 0.85]   # Light Gray
}
DEFAULT_COLOR = [1.00, 1.00, 1.00]

def render_alignment():
    alignment = AlignIO.read(INPUT_FASTA, "fasta")
    num_seqs = len(alignment)
    aln_len = alignment.get_alignment_length()

    start = max(0, START_POS - 1)
    end = min(aln_len, END_POS)
    sub_len = end - start

    res_per_block = int(np.ceil(sub_len / MAX_BLOCKS))
    num_blocks = int(np.ceil(sub_len / res_per_block))

    # Nature double-column width (183 mm = 7.205 in)
    fig_width_in = 183 / 25.4
    
    # Scaled height without extra space for a legend
    block_height_in = (num_seqs * 0.08 + 0.30) * (72 / 25.4) / 72
    fig_height_in = max(2.2, min(7.0, block_height_in * num_blocks + 0.3))

    fig, axes = plt.subplots(num_blocks, 1, figsize=(fig_width_in, fig_height_in), squeeze=False)

    for block_idx in range(num_blocks):
        ax = axes[block_idx, 0]
        b_start = start + block_idx * res_per_block
        b_end = min(start + (block_idx + 1) * res_per_block, end)
        cur_len = b_end - b_start

        # Build RGB matrix
        grid_rgb = np.zeros((num_seqs, cur_len, 3))
        for seq_idx, record in enumerate(alignment):
            sub_seq = str(record.seq[b_start:b_end]).upper()
            for x, char in enumerate(sub_seq):
                grid_rgb[seq_idx, x] = AA_PALETTE.get(char, DEFAULT_COLOR)

        # Render matrix
        ax.imshow(
            grid_rgb,
            aspect="auto",
            interpolation="nearest",
            extent=[-0.5, cur_len - 0.5, num_seqs - 0.5, -0.5]
        )

        # Taxa names
        ax.set_yticks(range(num_seqs))
        ax.set_yticklabels([r.id for r in alignment], fontsize=FONT_SIZE_LABEL, fontweight="normal")
        ax.tick_params(axis="y", length=0, pad=3, labelsize=FONT_SIZE_LABEL, labelcolor="#222222")

        # Top residue numbering
        tick_step = 20 if cur_len > 60 else 10
        x_ticks = np.arange(0, cur_len, tick_step)
        x_labels = [str(b_start + x + 1) for x in x_ticks]

        ax.set_xticks(x_ticks)
        ax.set_xticklabels(x_labels, fontsize=FONT_SIZE_NUM, color="#444444")
        ax.xaxis.tick_top()
        ax.tick_params(axis="x", length=2, pad=1, direction="out", color="#888888")

        # Clean frame
        for spine in ax.spines.values():
            spine.set_edgecolor("#CCCCCC")
            spine.set_linewidth(0.5)

        ax.set_xlim(-0.5, cur_len - 0.5)

    # Tight layout without bottom margin reservation for legend
    plt.subplots_adjust(left=0.18, right=0.96, top=0.94, bottom=0.04, hspace=0.50)
    plt.savefig(OUTPUT_PDF, format="pdf", bbox_inches="tight")
    print(f"Publication-ready AA alignment saved to: {OUTPUT_PDF}")

if __name__ == "__main__":
    render_alignment()