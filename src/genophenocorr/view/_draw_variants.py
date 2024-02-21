from itertools import cycle

import numpy as np
import matplotlib.pyplot as plt

from genophenocorr.model import Cohort, ProteinMetadata, TranscriptCoordinates, VariantEffect


#  BASIC DRAWING METHODS
def draw_rectangle(start_x, start_y, end_x, end_y, line_color='black', fill_color=None, line_width=1.0):
    rect = plt.Rectangle((start_x, start_y), end_x - start_x, end_y - start_y, edgecolor=line_color,
                         fill=fill_color is not None, linewidth=line_width, facecolor=fill_color)
    plt.gca().add_patch(rect)


def draw_line(start_x, start_y, end_x, end_y, line_color='black', line_width=1.0):
    plt.plot([start_x, end_x], [start_y, end_y], color=line_color, linewidth=line_width)


def draw_circle(center_x, center_y, radius, line_color='black', fill_color=None, line_width=1.0):
    circle = plt.Circle((center_x, center_y), radius, edgecolor=line_color, fill=fill_color is not None,
                        linewidth=line_width, facecolor=fill_color)
    plt.gca().add_patch(circle)


def draw_string(text, x, y, ha, va, color='black', fontsize=12, rotation=0):
    plt.text(x, y, text, fontsize=fontsize, color=color, ha=ha, va=va, rotation=rotation)


class VariantsVisualizer:
    def __init__(self):
        self.protein_track_color = '#d3d3d3'
        self.marker_colors = {  # TODO: @ielis change colors for variant effects as desired
            VariantEffect.TRANSCRIPT_ABLATION: "#ff0000",
            VariantEffect.SPLICE_ACCEPTOR_VARIANT: "#00ff00",
            VariantEffect.SPLICE_DONOR_VARIANT: "#0000ff",
            VariantEffect.STOP_GAINED: "#ff00ff",
            VariantEffect.FRAMESHIFT_VARIANT: "#ffff00",
            VariantEffect.STOP_LOST: "#00ffff",
            VariantEffect.START_LOST: "#ff9900",
            VariantEffect.TRANSCRIPT_AMPLIFICATION: "#9900ff",
            VariantEffect.INFRAME_INSERTION: "#ff0099",
            VariantEffect.INFRAME_DELETION: "#99ff00",
            VariantEffect.MISSENSE_VARIANT: "#00ff99",
            VariantEffect.PROTEIN_ALTERING_VARIANT: "#990000",
            VariantEffect.SPLICE_REGION_VARIANT: "#009900",
            VariantEffect.SPLICE_DONOR_5TH_BASE_VARIANT: "#009999",
            VariantEffect.SPLICE_DONOR_REGION_VARIANT: "#990099",
            VariantEffect.SPLICE_POLYPYRIMIDINE_TRACT_VARIANT: "#999900",
            VariantEffect.INCOMPLETE_TERMINAL_CODON_VARIANT: "#999999",
            VariantEffect.START_RETAINED_VARIANT: "#ffcc00",
            VariantEffect.STOP_RETAINED_VARIANT: "#ccff00",
            VariantEffect.SYNONYMOUS_VARIANT: "#00ccff",
            VariantEffect.CODING_SEQUENCE_VARIANT: "#ff00cc",
            VariantEffect.MATURE_MIRNA_VARIANT: "#cc00ff",
            VariantEffect.FIVE_PRIME_UTR_VARIANT: "#ff6600",
            VariantEffect.THREE_PRIME_UTR_VARIANT: "#6600ff",
            VariantEffect.NON_CODING_TRANSCRIPT_EXON_VARIANT: "#ff3366",
            VariantEffect.INTRON_VARIANT: "#3366ff",
            VariantEffect.NMD_TRANSCRIPT_VARIANT: "#ffcc99",
            VariantEffect.NON_CODING_TRANSCRIPT_VARIANT: "#cc99ff",
            VariantEffect.UPSTREAM_GENE_VARIANT: "#ff6633",
            VariantEffect.DOWNSTREAM_GENE_VARIANT: "#6633ff",
            VariantEffect.TFBS_ABLATION: "#cc3300",
            VariantEffect.TFBS_AMPLIFICATION: "#ccff66",
            VariantEffect.TF_BINDING_SITE_VARIANT: "#66ccff",
            VariantEffect.REGULATORY_REGION_ABLATION: "#ff3366",
            VariantEffect.REGULATORY_REGION_AMPLIFICATION: "#3366ff",
            VariantEffect.FEATURE_ELONGATION: "#ffcc33",
            VariantEffect.REGULATORY_REGION_VARIANT: "#ccff33",
            VariantEffect.FEATURE_TRUNCATION: "#33ccff",
            VariantEffect.INTERGENIC_VARIANT: "#ff0033",
            VariantEffect.SEQUENCE_VARIANT: "#33ff00",
        }
        self.feature_colors = ['red', 'orange', 'yellow']
        self.feature_outline_color = 'black'
        self.exon_colors = cycle(['blue', 'lightblue'])
        self.exon_outline_color = 'black'
        self.axis_color = 'black'

    def _draw_marker(self, x, min_y, max_y, circle_radius, color):
        draw_line(x, min_y, x, max_y, line_color=self.protein_track_color, line_width=0.5)
        draw_circle(x, max_y, circle_radius, line_color=self.protein_track_color, fill_color=color, line_width=0.5)

    def _marker_dim(self, marker_count, protein_track_y_max, marker_length=0.02, marker_radius=0.0025):
        radius = marker_radius + np.sqrt(marker_count - 1) * marker_radius
        length = protein_track_y_max + marker_length + np.sqrt(marker_count - 1) * marker_length
        return radius, length

    def _get_tx_anns(self, variants, tx_id):
        tx_anns = []
        for i, v in enumerate(variants):
            tx_ann = None
            for ann in v.tx_annotations:
                if ann.transcript_id == tx_id:
                    tx_ann = ann
                    break
            if tx_ann is None:
                raise ValueError(f'The transcript annotation for {tx_id} was not found!')
            else:
                tx_anns.append(tx_ann)

        return tx_anns

    def draw_fig(self, tx_coordinates: TranscriptCoordinates, protein_meta: ProteinMetadata, cohort: Cohort):
        tx_id = tx_coordinates.identifier
        protein_id = protein_meta.protein_id
        variants = cohort.all_variants
        tx_anns = self._get_tx_anns(variants, tx_id)

        exon_limits = np.array([(cds.start, cds.end) for cds in tx_coordinates.get_cds_regions()])
        feature_limits = np.array([(feature.info.start, feature.info.end) for feature in protein_meta.protein_features])
        variant_locations = np.array([[
            ann.protein_effect_location.start,
            ann.protein_effect_location.end]
            for ann in tx_anns
        ])
        variant_effects = np.array([(ann.variant_effects[0]) for ann in tx_anns])
        # count marker occurrences and remove duplicates
        variant_locations_counted_absolute, marker_counts = np.unique(variant_locations, axis=0, return_counts=True)
        variant_effect_colors = []
        for vl in variant_locations_counted_absolute:
            i = np.where(variant_locations == vl)[0][0]  # find index of unique variant loc in all locs to find effect
            effect = variant_effects[i]
            variant_effect_colors.append(self.marker_colors[effect])
        exon_labels = [f'{i + 1}' for i in range(len(exon_limits))]

        protein_track_x_min, protein_track_x_max = 0.15, 0.85
        protein_track_y_min, protein_track_y_max = 0.492, 0.508
        font_size = 12
        text_padding = 0.004

        plt.figure(figsize=(20, 20))

        min_x_absolute = min(np.min(feature_limits), np.min(exon_limits), np.min(variant_locations))
        max_x_absolute = max(np.max(feature_limits), np.max(exon_limits), np.max(variant_locations))

        max_marker_count = np.max(marker_counts)

        # normalize into [0, 1], leaving some space on the sides
        def preprocess(x_absolute):
            shifted_to_0_1 = ((x_absolute - min_x_absolute) / (max_x_absolute - min_x_absolute))
            relative_scale = (protein_track_x_max - protein_track_x_min)
            return shifted_to_0_1 * relative_scale + protein_track_x_min

        exon_limits_relative = preprocess(exon_limits)
        feature_limits_relative = preprocess(feature_limits)
        variant_locations_relative = preprocess(variant_locations_counted_absolute)

        # draw the protein track
        draw_rectangle(protein_track_x_min, protein_track_y_min, protein_track_x_max, protein_track_y_max,
                       line_color=self.protein_track_color, fill_color=self.protein_track_color, line_width=2.0)
        # x_axis
        x_axis_y = protein_track_y_min - 0.02
        x_axis_min_x, x_axis_max_x = protein_track_x_min, protein_track_x_max
        big_tick_length, small_tick_length = 0.01, 0.005
        draw_line(x_axis_min_x, x_axis_y, x_axis_max_x, x_axis_y, line_color=self.axis_color,
                  line_width=1.0)  # main line
        draw_line(x_axis_min_x, x_axis_y - big_tick_length, x_axis_min_x, x_axis_y, line_color=self.axis_color,
                  line_width=1.0)  # minimum tick
        draw_string(str(min_x_absolute), x_axis_min_x, x_axis_y - big_tick_length - text_padding, fontsize=font_size,
                    ha='center', va='top')
        draw_line(x_axis_max_x, x_axis_y - big_tick_length, x_axis_max_x, x_axis_y, line_color=self.axis_color,
                  line_width=1.0)  # max tick
        draw_string(str(max_x_absolute), x_axis_max_x, x_axis_y - big_tick_length - text_padding, fontsize=font_size,
                    ha='center', va='top')

        # y_axis
        y_axis_x = protein_track_x_min - 0.02
        y_axis_min_y = protein_track_y_max + 0.01
        _, y_axis_max_y = self._marker_dim(max_marker_count, protein_track_y_max)
        draw_line(y_axis_x, y_axis_min_y, y_axis_x, y_axis_max_y, line_color=self.axis_color, line_width=1.0)
        draw_line(y_axis_x - small_tick_length, y_axis_min_y, y_axis_x, y_axis_min_y, line_color=self.axis_color,
                  line_width=1.0)  # 0 tick
        draw_string("0", y_axis_x - small_tick_length - text_padding, y_axis_min_y, fontsize=font_size, ha='right',
                    va='center')
        draw_line(y_axis_x - small_tick_length, y_axis_max_y, y_axis_x, y_axis_max_y, line_color=self.axis_color,
                  line_width=1.0)  # max tick
        draw_string(str(max_marker_count), y_axis_x - small_tick_length - text_padding, y_axis_max_y,
                    fontsize=font_size, ha='right', va='center')
        draw_string("# Markers", y_axis_x - 0.05, (y_axis_min_y + y_axis_max_y) / 2, fontsize=font_size, ha='center',
                    va='center', rotation=90)  # x axis label

        # draw the exons (transcript track)
        transcript_y_min, transcript_y_max = 0.39, 0.43
        # iterate over pairs
        for exon_x_min, exon_x_max in exon_limits_relative:
            cur_color = next(self.exon_colors)
            draw_rectangle(exon_x_min, transcript_y_min, exon_x_max, transcript_y_max,
                           line_color=self.exon_outline_color, fill_color=cur_color, line_width=1.0)

        # draw variants
        marker_y_min = protein_track_y_max
        for marker in variant_locations_relative:
            marker_count = marker_counts[np.where(variant_locations_relative == marker)[0][0]]
            cur_radius, cur_length = self._marker_dim(marker_count, protein_track_y_max)
            self._draw_marker(marker, marker_y_min, cur_length, cur_radius, np.random.choice(self.marker_colors))

        # draw the features (protein track)
        feature_y_min, feature_y_max = 0.485, 0.515
        for feature_x_min, feature_x_max in feature_limits_relative:
            draw_rectangle(feature_x_min, feature_y_min, feature_x_max, feature_y_max,
                           line_color=self.feature_outline_color,
                           fill_color=np.random.choice(self.feature_colors), line_width=1.0)

        plt.xlim(0, 1)
        plt.ylim(0.3, 0.7)
        plt.gca().set_aspect('equal')
        plt.axis('off')
        plt.title(f'[Working title:] transcript: {tx_id}, '
                  f'protein: {protein_id},'
                  f'protein name: {protein_meta.label}')
        plt.show()
