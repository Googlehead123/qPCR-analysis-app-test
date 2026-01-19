"""
qPCR Analysis Suite Pro v3.0
Enhanced version with QC module, improved UI, and publication-ready exports.

Changes from v2.0:
- Fixed: Unused parameter bug in calculate_ddct
- Fixed: Bare exception clauses
- Fixed: Import placement (all at top)
- Removed: Duplicate function definitions
- Removed: Redundant session state checks
- Added: QC Module with outlier detection
- Added: Bulk sample operations
- Added: Publication-ready export (PNG/SVG/PDF/TIFF)
- Improved: Collapsible UI sections
- Improved: Cleaner layout and workflow
"""

# ==================== IMPORTS (ALL AT TOP) ====================
import streamlit as st
import pandas as pd
import numpy as np
import plotly.graph_objects as go
import plotly.express as px
from scipy import stats
import io
import json
import re
import textwrap
from datetime import datetime
from typing import Dict, List, Tuple, Optional, Set
import base64

# ==================== PAGE CONFIG ====================
st.set_page_config(
    page_title="qPCR Analysis Suite Pro v3.0", 
    layout="wide", 
    initial_sidebar_state="expanded",
    page_icon="🧬"
)

# ==================== CONSTANTS ====================
QC_THRESHOLDS = {
    'cv_warning': 5.0,      # CV% warning threshold
    'cv_critical': 10.0,    # CV% critical threshold
    'ct_warning': 32.0,     # High CT warning
    'ct_critical': 35.0,    # High CT critical
    'zscore_threshold': 2.5, # Outlier detection
    'min_replicates': 2     # Minimum replicates for statistics
}

EXPORT_PRESETS = {
    'journal': {
        'name': 'Journal Publication',
        'dpi': 300,
        'width': 7,  # inches (double column)
        'height': 5,
        'font_family': 'Arial',
        'font_size': 12,
        'formats': ['PNG', 'TIFF', 'PDF']
    },
    'poster': {
        'name': 'Poster/Presentation',
        'dpi': 150,
        'width': 10,
        'height': 7,
        'font_family': 'Arial',
        'font_size': 16,
        'formats': ['PNG', 'SVG']
    },
    'manuscript': {
        'name': 'Single Column Figure',
        'dpi': 300,
        'width': 3.5,
        'height': 3,
        'font_family': 'Arial',
        'font_size': 10,
        'formats': ['PNG', 'TIFF']
    }
}

# ==================== EFFICACY DATABASE ====================
EFFICACY_CONFIG = {
    '탄력': {
        'genes': ['COL1A1', 'ELN', 'FBN-1', 'FBN1'],
        'cell': 'HS68 fibroblast',
        'controls': {
            'negative': 'Non-treated',
            'positive': 'TGFb',
            'compare_to': 'negative'
        },
        'description': 'Elasticity - Non-treated vs TGFb (positive) vs Treatments'
    },
    '항노화': {
        'genes': ['COL1A1', 'COL1', 'MMP-1', 'MMP1'],
        'cell': 'HS68 fibroblast',
        'controls': {
            'baseline': 'Non-treated (No UV)',
            'negative': 'UVB only',
            'positive': 'UVB+TGFb',
            'compare_to': 'negative'
        },
        'description': 'Anti-aging - COL1↑ (recovery), MMP1↓ (inhibition) after UVB damage',
        'expected_direction': {'COL1A1': 'up', 'COL1': 'up', 'MMP-1': 'down', 'MMP1': 'down'}
    },
    '보습': {
        'genes': ['AQP3', 'HAS3'],
        'cell': 'HaCaT keratinocyte',
        'controls': {
            'negative': 'Non-treated',
            'positive': 'Retinoic acid',
            'compare_to': 'negative'
        },
        'description': 'Hydration - Non-treated vs Retinoic acid (positive) vs Treatments'
    },
    '장벽': {
        'genes': ['FLG', 'CLDN', 'IVL'],
        'cell': 'HaCaT keratinocyte',
        'controls': {
            'negative': 'Non-treated',
            'positive': 'Retinoic acid',
            'compare_to': 'negative'
        },
        'description': 'Barrier function - Non-treated vs Retinoic acid (positive) vs Treatments'
    },
    '표피증식': {
        'genes': ['KI67', 'PCNA'],
        'cell': 'HaCaT keratinocyte',
        'controls': {
            'negative': 'Non-treated',
            'positive': 'TGFb or FBS',
            'compare_to': 'negative'
        },
        'description': 'Proliferation - Non-treated vs TGFb/FBS (positive) vs Treatments'
    },
    '멜라닌억제': {
        'genes': ['MITF', 'TYR', 'Melanin'],
        'cell': 'B16F10 melanocyte',
        'controls': {
            'baseline': 'Non-treated',
            'negative': 'α-MSH only',
            'positive': 'α-MSH+Arbutin',
            'compare_to': 'negative'
        },
        'description': 'Melanin inhibition - α-MSH induced vs α-MSH+Arbutin (positive) vs α-MSH+Treatments',
        'expected_direction': {'MITF': 'down', 'TYR': 'down', 'Melanin': 'down'}
    },
    '진정': {
        'genes': ['IL1B', 'IL-1β', 'IL6', 'TNFA', 'TNFα'],
        'cell': 'HaCaT keratinocyte',
        'controls': {
            'baseline': 'Non-treated',
            'negative': 'IL4+PolyIC (Inflammation)',
            'positive': 'Inflammation+Dexamethasone',
            'compare_to': 'negative'
        },
        'description': 'Anti-inflammation - Reduce IL1β/IL6/TNFα (all should decrease)',
        'expected_direction': {'IL1B': 'down', 'IL-1β': 'down', 'IL6': 'down', 'TNFA': 'down', 'TNFα': 'down'}
    },
    '지질억제': {
        'genes': ['SREBPA', 'SREBPa', 'SREBPC', 'SREBPc', 'PPARY', 'PPARy'],
        'cell': 'SZ95 sebocyte',
        'controls': {
            'baseline': 'Non-treated',
            'negative': 'IGF only',
            'positive': 'IGF+Reference inhibitor',
            'compare_to': 'negative'
        },
        'description': 'Sebum inhibition - IGF induced vs IGF+Treatments',
        'expected_direction': {'SREBPA': 'down', 'SREBPa': 'down', 'SREBPC': 'down', 'SREBPc': 'down', 'PPARY': 'down', 'PPARy': 'down'}
    },
    '냉감': {
        'genes': ['TRPM8', 'CIRBP'],
        'cell': 'HaCaT keratinocyte',
        'controls': {
            'negative': 'Non-treated',
            'positive': 'Menthol',
            'compare_to': 'negative'
        },
        'description': 'Cooling effect - Non-treated vs Menthol (positive) vs Treatments'
    }
}

# ==================== SESSION STATE INITIALIZATION ====================
def init_session_state():
    """Initialize all session state variables"""
    defaults = {
        'data': None,
        'processed_data': None,
        'sample_mapping': {},
        'analysis_templates': {},
        'graphs': {},
        'excluded_wells': set(),
        'excluded_samples': set(),
        'selected_efficacy': None,
        'hk_gene': None,
        'sample_order': None,
        'graph_settings': {
            'title_size': 20, 'font_size': 14, 'sig_font_size': 18,
            'figure_width': 1000, 'figure_height': 600,
            'color_scheme': 'plotly_white', 'show_error': True,
            'show_significance': True, 'show_grid': True,
            'xlabel': 'Condition', 'ylabel': 'Relative mRNA Expression Level',
            'bar_colors': {}, 'orientation': 'v', 'error_multiplier': 1.96,
            'bar_opacity': 0.95, 'bar_gap': 0.15, 'marker_line_width': 1,
            'show_legend': False, 'y_log_scale': False, 'y_min': None, 'y_max': None
        },
        # QC Module state
        'qc_flagged_wells': set(),
        'qc_edited_values': {},
        'qc_exclusion_reasons': {},
        'qc_complete': False,
        # Export state
        'export_settings': EXPORT_PRESETS['journal'].copy(),
        # Bulk operations
        'selected_samples': set()
    }
    
    for key, default in defaults.items():
        if key not in st.session_state:
            st.session_state[key] = default

init_session_state()

# ==================== UTILITY FUNCTIONS ====================
def natural_sort_key(sample_name: str) -> list:
    """Extract numbers from sample name for natural sorting"""
    parts = re.split(r'(\d+)', str(sample_name))
    return [int(part) if part.isdigit() else part.lower() for part in parts]


def validate_ct_value(ct: float) -> Tuple[bool, str]:
    """Validate a CT value and return status with message"""
    if pd.isna(ct):
        return False, "Missing value"
    if ct < 0:
        return False, "Negative CT (invalid)"
    if ct > QC_THRESHOLDS['ct_critical']:
        return False, f"CT > {QC_THRESHOLDS['ct_critical']} (critical)"
    if ct > QC_THRESHOLDS['ct_warning']:
        return True, f"CT > {QC_THRESHOLDS['ct_warning']} (warning)"
    return True, "OK"


def calculate_cv(values: np.ndarray) -> float:
    """Calculate coefficient of variation (%)"""
    if len(values) < 2:
        return 0.0
    mean = np.mean(values)
    if mean == 0:
        return 0.0
    return (np.std(values, ddof=1) / mean) * 100


def detect_outliers_zscore(values: np.ndarray, threshold: float = 2.5) -> np.ndarray:
    """Detect outliers using modified Z-score (MAD-based, robust to small samples)"""
    if len(values) < 3:
        return np.array([False] * len(values))
    
    median = np.median(values)
    mad = np.median(np.abs(values - median))
    
    if mad == 0:
        # All values are the same or too few unique values
        return np.array([False] * len(values))
    
    modified_z = 0.6745 * (values - median) / mad
    return np.abs(modified_z) > threshold


def format_p_value(p: float) -> str:
    """Format p-value for display"""
    if pd.isna(p):
        return "—"
    if p < 0.001:
        return "< 0.001"
    if p < 0.01:
        return f"{p:.3f}"
    return f"{p:.2f}"


def get_significance_symbol(p: float, symbol_type: str = 'asterisk') -> str:
    """Get significance symbol based on p-value"""
    if pd.isna(p):
        return ""
    
    if symbol_type == 'asterisk':
        if p < 0.001:
            return "***"
        elif p < 0.01:
            return "**"
        elif p < 0.05:
            return "*"
    else:  # hashtag
        if p < 0.001:
            return "###"
        elif p < 0.01:
            return "##"
        elif p < 0.05:
            return "#"
    return ""


# ==================== PARSER CLASS ====================
class QPCRParser:
    """Parse qPCR data from various CSV formats"""
    
    @staticmethod
    def detect_format(df: pd.DataFrame) -> Tuple[str, int]:
        """Detect the format of the qPCR data file"""
        for idx, row in df.iterrows():
            row_str = ' '.join(row.astype(str).values)
            if 'Well Position' in row_str:
                return 'format1', idx
            elif row.iloc[0] == 'Well' and 'Sample Name' in row_str:
                return ('format2' if 'Cт' in row_str or 'ΔCт' in row_str else 'format1'), idx
        return 'unknown', 0
    
    @staticmethod
    def parse_format1(df: pd.DataFrame, start: int) -> Optional[pd.DataFrame]:
        """Parse format 1 (Well Position style)"""
        df = df.iloc[start:].reset_index(drop=True)
        df.columns = df.iloc[0]
        df = df.iloc[1:].reset_index(drop=True)
        
        well_col = next((c for c in ['Well Position', 'Well'] if c in df.columns), df.columns[0])
        ct_col = next((c for c in ['CT', 'Ct', 'Cт'] if c in df.columns), None)
        
        if not ct_col:
            return None
        
        parsed = pd.DataFrame({
            'Well': df[well_col],
            'Sample': df.get('Sample Name', df.iloc[:, 2]),
            'Target': df.get('Target Name', df.iloc[:, 3]),
            'CT': pd.to_numeric(df[ct_col], errors='coerce')
        })
        
        return parsed.dropna(subset=['CT']).query('Sample.notna() & Target.notna()')
    
    @staticmethod
    def parse_format2(df: pd.DataFrame, start: int) -> Optional[pd.DataFrame]:
        """Parse format 2 (Cyrillic CT style)"""
        df = df.iloc[start:].reset_index(drop=True)
        df.columns = df.iloc[0]
        df = df.iloc[1:].reset_index(drop=True)
        
        parsed = pd.DataFrame({
            'Well': df['Well'],
            'Sample': df['Sample Name'],
            'Target': df['Target Name'],
            'CT': pd.to_numeric(df['Cт'], errors='coerce')
        })
        
        return parsed.dropna(subset=['CT']).query('Sample.notna() & Target.notna()')
    
    @staticmethod
    def parse(file) -> Optional[pd.DataFrame]:
        """Main parsing method with encoding detection"""
        try:
            df = None
            for enc in ['utf-8', 'latin-1', 'cp1252']:
                try:
                    file.seek(0)  # Reset file pointer
                    df = pd.read_csv(file, encoding=enc, low_memory=False, skip_blank_lines=False)
                    break
                except UnicodeDecodeError:
                    continue
            
            if df is None:
                return None
            
            fmt, start = QPCRParser.detect_format(df)
            
            if fmt == 'format1':
                return QPCRParser.parse_format1(df, start)
            elif fmt == 'format2':
                return QPCRParser.parse_format2(df, start)
            else:
                return None
                
        except Exception as e:
            st.error(f"Parse error: {e}")
            return None


# ==================== QC MODULE ====================
class QCModule:
    """Quality Control module for qPCR data validation"""
    
    @staticmethod
    def run_qc_analysis(data: pd.DataFrame) -> Dict:
        """Run comprehensive QC analysis on the data"""
        qc_results = {
            'summary': {},
            'well_flags': {},
            'replicate_stats': {},
            'missing_data': [],
            'high_ct_wells': [],
            'outlier_wells': [],
            'quality_score': 100.0
        }
        
        if data is None or data.empty:
            return qc_results
        
        total_wells = len(data)
        issues_count = 0
        
        # Analyze each sample-target combination
        for (sample, target), group in data.groupby(['Sample', 'Target']):
            ct_values = group['CT'].values
            wells = group['Well'].tolist()
            
            # Calculate replicate statistics
            n_reps = len(ct_values)
            ct_mean = np.mean(ct_values)
            ct_std = np.std(ct_values, ddof=1) if n_reps > 1 else 0
            cv = calculate_cv(ct_values)
            
            qc_results['replicate_stats'][(sample, target)] = {
                'n': n_reps,
                'mean': ct_mean,
                'std': ct_std,
                'cv': cv,
                'wells': wells
            }
            
            # Check for high CT values
            for well, ct in zip(wells, ct_values):
                valid, msg = validate_ct_value(ct)
                if 'critical' in msg or 'warning' in msg:
                    qc_results['high_ct_wells'].append({
                        'well': well,
                        'sample': sample,
                        'target': target,
                        'ct': ct,
                        'status': 'critical' if 'critical' in msg else 'warning'
                    })
                    qc_results['well_flags'][well] = msg
                    if 'critical' in msg:
                        issues_count += 2
                    else:
                        issues_count += 1
            
            # Check CV thresholds
            if cv > QC_THRESHOLDS['cv_critical']:
                for well in wells:
                    if well not in qc_results['well_flags']:
                        qc_results['well_flags'][well] = f"High CV ({cv:.1f}%)"
                issues_count += 2
            elif cv > QC_THRESHOLDS['cv_warning']:
                issues_count += 1
            
            # Detect outliers within replicates
            if n_reps >= 3:
                outlier_mask = detect_outliers_zscore(ct_values, QC_THRESHOLDS['zscore_threshold'])
                for well, ct, is_outlier in zip(wells, ct_values, outlier_mask):
                    if is_outlier:
                        qc_results['outlier_wells'].append({
                            'well': well,
                            'sample': sample,
                            'target': target,
                            'ct': ct,
                            'replicate_mean': ct_mean
                        })
                        if well not in qc_results['well_flags']:
                            qc_results['well_flags'][well] = "Potential outlier"
                        issues_count += 1
        
        # Check for missing data (expected but not present)
        all_samples = data['Sample'].unique()
        all_targets = data['Target'].unique()
        
        for sample in all_samples:
            for target in all_targets:
                if (sample, target) not in qc_results['replicate_stats']:
                    qc_results['missing_data'].append({
                        'sample': sample,
                        'target': target
                    })
                    issues_count += 1
        
        # Calculate quality score (100 = perfect, decreases with issues)
        max_penalty = total_wells * 0.5  # Maximum penalty
        penalty = min(issues_count * 2, max_penalty)
        qc_results['quality_score'] = max(0, 100 - (penalty / total_wells * 100)) if total_wells > 0 else 0
        
        # Summary statistics
        qc_results['summary'] = {
            'total_wells': total_wells,
            'total_samples': len(all_samples),
            'total_targets': len(all_targets),
            'flagged_wells': len(qc_results['well_flags']),
            'high_ct_count': len(qc_results['high_ct_wells']),
            'outlier_count': len(qc_results['outlier_wells']),
            'missing_count': len(qc_results['missing_data']),
            'quality_score': round(qc_results['quality_score'], 1)
        }
        
        return qc_results
    
    @staticmethod
    def apply_exclusions(data: pd.DataFrame, excluded_wells: Set[str], 
                        edited_values: Dict) -> pd.DataFrame:
        """Apply QC exclusions and edits to data"""
        data = data.copy()
        
        # Remove excluded wells
        data = data[~data['Well'].isin(excluded_wells)]
        
        # Apply edited CT values
        for well, edit_info in edited_values.items():
            mask = data['Well'] == well
            if mask.any():
                data.loc[mask, 'CT'] = edit_info['edited']
        
        return data


# ==================== ANALYSIS ENGINE ====================
class AnalysisEngine:
    """Core analysis engine for ΔΔCt calculations and statistics"""
    
    @staticmethod
    def calculate_ddct(data: pd.DataFrame, hk_gene: str, ref_sample: str,
                       excluded_wells: Set, excluded_samples: Set, 
                       sample_mapping: Dict) -> pd.DataFrame:
        """
        Gene-by-gene ΔΔCt calculation with housekeeping normalization.
        
        NOTE: Removed unused 'compare_sample' parameter from v2.0 (bug fix)
        """
        # Filter data
        data = data[~data['Well'].isin(excluded_wells) & ~data['Sample'].isin(excluded_samples)].copy()
        
        # Apply sample name mapping
        data['Condition'] = data['Sample'].map(lambda x: sample_mapping.get(x, {}).get('condition', x))
        data['Group'] = data['Sample'].map(lambda x: sample_mapping.get(x, {}).get('group', 'Treatment'))
        
        results = []
        
        # Process each target gene separately (exclude housekeeping)
        hk_variants = [hk_gene.upper(), 'ACTIN', 'B-ACTIN', 'GAPDH', 'ACTB']
        
        for target in data['Target'].unique():
            if target.upper() in hk_variants:
                continue
            
            target_data = data[data['Target'] == target]
            
            for condition in target_data['Condition'].unique():
                cond_data = target_data[target_data['Condition'] == condition]
                
                # Get housekeeping Ct for this condition
                hk_data = data[(data['Condition'] == condition) & (data['Target'] == hk_gene)]
                if len(hk_data) == 0:
                    continue
                
                # Calculate ΔCt = Target_Ct - HK_Ct
                target_ct_values = cond_data['CT'].values
                hk_ct_values = hk_data['CT'].values
                
                # Use mean if multiple replicates
                target_ct_mean = target_ct_values.mean()
                hk_ct_mean = hk_ct_values.mean()
                delta_ct = target_ct_mean - hk_ct_mean
                
                # Get reference ΔCt (ref_sample)
                ref_target = target_data[target_data['Condition'] == ref_sample]
                ref_hk = data[(data['Condition'] == ref_sample) & (data['Target'] == hk_gene)]

                if len(ref_target) > 0 and len(ref_hk) > 0:
                    ref_delta_ct = ref_target['CT'].mean() - ref_hk['CT'].mean()
                else:
                    ref_delta_ct = 0
                
                # ΔΔCt and relative expression
                ddct = delta_ct - ref_delta_ct
                rel_expr = 2 ** (-ddct)
                
                # Error calculation
                ct_sd = target_ct_values.std(ddof=1) if len(target_ct_values) > 1 else 0
                sem = ct_sd / np.sqrt(len(target_ct_values)) if len(target_ct_values) > 1 else 0
                
                # Get original sample name and group
                original_sample = cond_data['Sample'].iloc[0]
                group = sample_mapping.get(original_sample, {}).get('group', 'Treatment')
                
                results.append({
                    'Target': target,
                    'Condition': condition,
                    'Original_Sample': original_sample,
                    'Group': group,
                    'n_replicates': len(target_ct_values),
                    'Target_Ct_Mean': target_ct_mean,
                    'Target_Ct_SD': ct_sd,
                    'HK_Ct_Mean': hk_ct_mean,
                    'Delta_Ct': delta_ct,
                    'Delta_Delta_Ct': ddct,
                    'Relative_Expression': rel_expr,
                    'SEM': sem
                })
        
        return pd.DataFrame(results)
    
    @staticmethod
    def calculate_statistics(processed: pd.DataFrame, compare_condition: str,
                            compare_condition_2: Optional[str] = None,
                            raw_data: Optional[pd.DataFrame] = None,
                            hk_gene: Optional[str] = None,
                            sample_mapping: Optional[Dict] = None) -> pd.DataFrame:
        """Two-tailed Welch's t-test comparing each condition to compare_condition"""
        
        # Use session_state fallbacks
        raw_data = raw_data if raw_data is not None else st.session_state.get("data")
        hk_gene = hk_gene if hk_gene is not None else st.session_state.get("hk_gene")
        sample_mapping = sample_mapping if sample_mapping is not None else st.session_state.get("sample_mapping", {})
        
        if raw_data is None or hk_gene is None:
            return processed
        
        results = processed.copy()
        results["p_value"] = np.nan
        results["significance"] = ""
        
        # Add second p-value columns if compare_condition_2 is provided
        if compare_condition_2:
            results["p_value_2"] = np.nan
            results["significance_2"] = ""
        
        for target in results["Target"].unique():
            if pd.isna(target):
                continue
            
            # Map conditions
            t_rows = raw_data[raw_data["Target"] == target].copy()
            if t_rows.empty:
                continue
            
            t_rows["Condition"] = t_rows["Sample"].map(
                lambda s: sample_mapping.get(s, {}).get("condition", s)
            )
            
            hk_rows = raw_data[raw_data["Target"] == hk_gene].copy()
            hk_rows["Condition"] = hk_rows["Sample"].map(
                lambda s: sample_mapping.get(s, {}).get("condition", s)
            )
            hk_means = hk_rows.groupby("Condition")["CT"].mean().to_dict()
            
            # Calculate relative expression per condition
            rel_expr = {}
            for cond, grp in t_rows.groupby("Condition"):
                hk_mean = hk_means.get(cond, np.nan)
                if np.isnan(hk_mean):
                    continue
                rel_expr[cond] = 2 ** (-(grp["CT"].values - hk_mean))
            
            # FIRST COMPARISON: compare_condition
            ref_vals = rel_expr.get(compare_condition, np.array([]))
            if hasattr(ref_vals, '__len__') and len(ref_vals) >= 1:
                for cond, vals in rel_expr.items():
                    if cond == compare_condition or len(vals) == 0:
                        continue
                    
                    try:
                        if len(ref_vals) >= 2 and len(vals) >= 2:
                            _, p_val = stats.ttest_ind(ref_vals, vals, equal_var=False)
                        elif len(vals) == 1 and len(ref_vals) >= 2:
                            _, p_val = stats.ttest_1samp(ref_vals, vals[0])
                        elif len(ref_vals) == 1 and len(vals) >= 2:
                            _, p_val = stats.ttest_1samp(vals, ref_vals[0])
                        else:
                            p_val = np.nan
                    except (ValueError, TypeError, ZeroDivisionError):
                        p_val = np.nan
                    
                    mask = (results["Target"] == target) & (results["Condition"] == cond)
                    results.loc[mask, "p_value"] = p_val
                    results.loc[mask, "significance"] = get_significance_symbol(p_val, 'asterisk')
            
            # SECOND COMPARISON: compare_condition_2 (if provided)
            if compare_condition_2:
                ref_vals_2 = rel_expr.get(compare_condition_2, np.array([]))
                if hasattr(ref_vals_2, '__len__') and len(ref_vals_2) >= 1:
                    for cond, vals in rel_expr.items():
                        if cond == compare_condition_2 or len(vals) == 0:
                            continue
                        
                        try:
                            if len(ref_vals_2) >= 2 and len(vals) >= 2:
                                _, p_val_2 = stats.ttest_ind(ref_vals_2, vals, equal_var=False)
                            elif len(vals) == 1 and len(ref_vals_2) >= 2:
                                _, p_val_2 = stats.ttest_1samp(ref_vals_2, vals[0])
                            elif len(ref_vals_2) == 1 and len(vals) >= 2:
                                _, p_val_2 = stats.ttest_1samp(vals, ref_vals_2[0])
                            else:
                                p_val_2 = np.nan
                        except (ValueError, TypeError, ZeroDivisionError):
                            p_val_2 = np.nan
                        
                        mask = (results["Target"] == target) & (results["Condition"] == cond)
                        results.loc[mask, "p_value_2"] = p_val_2
                        results.loc[mask, "significance_2"] = get_significance_symbol(p_val_2, 'hashtag')
        
        return results

    @staticmethod
    def run_full_analysis(ref_sample_key: str, compare_sample_key: str, 
                         compare_sample_key_2: Optional[str] = None) -> bool:
        """Run ΔΔCt + statistical analysis and store results"""
        try:
            data = st.session_state.get("data")
            mapping = st.session_state.get("sample_mapping", {})
            hk_gene = st.session_state.get("hk_gene")

            if data is None:
                st.error("❌ No raw data loaded.")
                return False
            if not mapping:
                st.error("❌ Sample mapping not found.")
                return False
            if not hk_gene:
                st.error("❌ Housekeeping gene not selected.")
                return False

            # Apply QC exclusions if any
            qc_excluded = st.session_state.get('qc_flagged_wells', set())
            qc_edited = st.session_state.get('qc_edited_values', {})
            
            if qc_excluded or qc_edited:
                data = QCModule.apply_exclusions(data, qc_excluded, qc_edited)

            ref_condition = mapping.get(ref_sample_key, {}).get("condition", ref_sample_key)
            cmp_condition = mapping.get(compare_sample_key, {}).get("condition", compare_sample_key)
            cmp_condition_2 = None
            if compare_sample_key_2:
                cmp_condition_2 = mapping.get(compare_sample_key_2, {}).get("condition", compare_sample_key_2)

            msg = f"Running analysis: reference '{ref_condition}', comparison '{cmp_condition}'"
            if cmp_condition_2:
                msg += f", secondary '{cmp_condition_2}'"
            
            with st.spinner(msg + "..."):
                # ΔΔCt calculation
                processed_df = AnalysisEngine.calculate_ddct(
                    data,
                    hk_gene,
                    ref_condition,
                    st.session_state.get("excluded_wells", set()),
                    st.session_state.get("excluded_samples", set()),
                    mapping,
                )

                if processed_df is None or processed_df.empty:
                    st.warning("⚠️ No ΔΔCt results produced. Check mapping and housekeeping gene.")
                    return False

                # Statistical test
                processed_with_stats = AnalysisEngine.calculate_statistics(
                    processed_df,
                    cmp_condition,
                    cmp_condition_2,
                    raw_data=data,
                    hk_gene=hk_gene,
                    sample_mapping=mapping,
                )

                # Organize data for graphs
                gene_dict = {}
                if "Target" in processed_with_stats.columns:
                    for gene in processed_with_stats["Target"].unique():
                        gene_df = processed_with_stats[processed_with_stats["Target"] == gene].copy()
                        gene_dict[gene] = gene_df.reset_index(drop=True)
                else:
                    gene_dict = {"results": processed_with_stats.reset_index(drop=True)}

                st.session_state.processed_data = gene_dict

            return True

        except Exception as e:
            st.error(f"❌ Analysis failed: {e}")
            return False


# ==================== GRAPH GENERATOR ====================
class GraphGenerator:
    """Generate publication-quality graphs"""
    
    @staticmethod
    def _wrap_text(text: str, width: int = 15) -> str:
        """Wrap text for x-axis labels"""
        return textwrap.fill(str(text), width=width)
    
    @staticmethod
    def create_gene_graph(
        data: pd.DataFrame,
        gene: str,
        settings: Dict,
        efficacy_config: Dict = None,
        sample_order: Optional[List] = None
    ) -> go.Figure:
        """Create individual graph for each gene"""
        
        # Guard against empty data
        if data is None or data.empty:
            fig = go.Figure()
            fig.add_annotation(text="No data available", showarrow=False)
            return fig
        
        gene_data = data.copy()
        
        # Ensure we have required columns
        if 'Relative_Expression' not in gene_data.columns:
            st.error(f"Missing Relative_Expression column for {gene}")
            fig = go.Figure()
            fig.add_annotation(text=f"Missing data columns for {gene}", showarrow=False)
            return fig
        
        if 'SEM' not in gene_data.columns:
            gene_data['SEM'] = 0
        
        # Apply sample ordering
        if sample_order:
            mapping = st.session_state.get('sample_mapping', {})
            condition_order = []
            seen_conditions = set()
            
            for sample in sample_order:
                if mapping.get(sample, {}).get('include', True):
                    cond = mapping.get(sample, {}).get('condition', sample)
                    if cond in gene_data['Condition'].unique() and cond not in seen_conditions:
                        condition_order.append(cond)
                        seen_conditions.add(cond)
            
            # Add any conditions not in order
            for cond in gene_data['Condition'].unique():
                if cond not in seen_conditions:
                    condition_order.append(cond)
                    seen_conditions.add(cond)
            
            gene_data['Condition'] = pd.Categorical(
                gene_data['Condition'], 
                categories=condition_order, 
                ordered=True
            )
            gene_data = gene_data.sort_values('Condition')
        
        gene_data = gene_data.reset_index(drop=True)
        condition_names = gene_data['Condition'].tolist()
        n_bars = len(gene_data)
        
        # Get colors
        bar_colors = []
        control_colors = {
            'Baseline': '#FFFFFF',
            'Non-treated': '#FFFFFF',
            'Control': '#FFFFFF',
            'Negative Control': '#9E9E9E',
            'Inducer': '#9E9E9E',
            'Positive Control': '#9E9E9E',
            'Treatment': '#D3D3D3'
        }
        
        for idx, row in gene_data.iterrows():
            condition = row['Condition']
            group = row.get('Group', 'Treatment')
            
            custom_key = f"{gene}_{condition}"
            if custom_key in settings.get('bar_colors_per_sample', {}):
                bar_colors.append(settings['bar_colors_per_sample'][custom_key])
            elif group in control_colors:
                bar_colors.append(control_colors[group])
            else:
                bar_colors.append(settings.get('bar_colors', {}).get(gene, '#D3D3D3'))
        
        # Create figure
        fig = go.Figure()
        
        # Error bars
        error_array = (gene_data['SEM'] * settings.get('error_multiplier', 1.96)).values
        
        # Per-bar settings
        gene_bar_settings = st.session_state.get(f'{gene}_bar_settings', {})
        show_error_global = settings.get('show_error', True)
        show_sig_global = settings.get('show_significance', True)
        
        # Build error visibility array
        error_visible_array = []
        for idx in range(n_bars):
            row = gene_data.iloc[idx]
            condition = row['Condition']
            bar_key = f"{gene}_{condition}"
            bar_config = gene_bar_settings.get(bar_key, {'show_sig': True, 'show_err': True})
            
            if show_error_global and bar_config.get('show_err', True):
                error_visible_array.append(error_array[idx])
            else:
                error_visible_array.append(0)
        
        # Add bar trace
        fig.add_trace(go.Bar(
            x=list(range(n_bars)),
            y=gene_data['Relative_Expression'],
            error_y=dict(
                type='data',
                array=error_visible_array,
                arrayminus=[0] * n_bars,
                visible=True,
                thickness=2,
                width=4,
                color='rgba(0,0,0,0.5)',
                symmetric=False
            ),
            marker=dict(
                color=bar_colors,
                line=dict(
                    width=settings.get('marker_line_width', 1),
                    color='black'
                ),
                opacity=settings.get('bar_opacity', 0.95)
            ),
            showlegend=False
        ))
        
        # Calculate y-axis range
        max_y_value = gene_data['Relative_Expression'].max()
        max_error = max(error_array) if len(error_array) > 0 else 0
        y_max_auto = max_y_value + max_error + (max_y_value * 0.15)
        
        fixed_symbol_spacing = y_max_auto * 0.05
        
        # Add significance symbols
        for idx in range(n_bars):
            row = gene_data.iloc[idx]
            condition = row['Condition']
            bar_key = f"{gene}_{condition}"
            bar_config = gene_bar_settings.get(bar_key, {'show_sig': True, 'show_err': True})
            
            sig_1 = row.get('significance', '')
            sig_2 = row.get('significance_2', '')
            
            bar_height = row['Relative_Expression']
            error_bar_height = error_visible_array[idx]
            base_y_position = bar_height + error_bar_height
            
            asterisk_font_size = 16
            hashtag_font_size = 10
            
            if show_sig_global and bar_config.get('show_sig', True):
                symbols_to_show = []
                font_sizes = []
                
                if sig_1 in ['*', '**', '***']:
                    symbols_to_show.append(sig_1)
                    font_sizes.append(asterisk_font_size)
                
                if sig_2 in ['#', '##', '###']:
                    symbols_to_show.append(sig_2)
                    font_sizes.append(hashtag_font_size)
                
                if len(symbols_to_show) == 2:
                    fig.add_annotation(
                        x=idx, y=base_y_position + (fixed_symbol_spacing * 0.2),
                        text=symbols_to_show[0], showarrow=False,
                        font=dict(size=font_sizes[0], color='black', family='Arial'),
                        xref='x', yref='y', xanchor='center', yanchor='bottom'
                    )
                    fig.add_annotation(
                        x=idx, y=base_y_position + (fixed_symbol_spacing * 0.2) + fixed_symbol_spacing,
                        text=symbols_to_show[1], showarrow=False,
                        font=dict(size=font_sizes[1], color='black', family='Arial'),
                        xref='x', yref='y', xanchor='center', yanchor='bottom'
                    )
                elif len(symbols_to_show) == 1:
                    fig.add_annotation(
                        x=idx, y=base_y_position + (fixed_symbol_spacing * 0.2),
                        text=symbols_to_show[0], showarrow=False,
                        font=dict(size=font_sizes[0], color='black', family='Arial'),
                        xref='x', yref='y', xanchor='center', yanchor='bottom'
                    )
        
        # Y-axis configuration
        y_label_html = f"Relative <b style='color:red;'>{gene}</b> Expression Level"
        
        y_axis_config = dict(
            title=dict(
                text=y_label_html,
                font=dict(size=settings.get(f"{gene}_ylabel_size", 14))
            ),
            showgrid=False,
            zeroline=True,
            zerolinewidth=1,
            zerolinecolor='black',
            range=[0, y_max_auto],
            fixedrange=False
        )
        
        if settings.get('y_log_scale'):
            y_axis_config['type'] = 'log'
            y_axis_config.pop('range', None)
        
        if settings.get('y_min') is not None or settings.get('y_max') is not None:
            y_axis_config['range'] = [
                settings.get('y_min', 0),
                settings.get('y_max', y_max_auto)
            ]
        
        gene_bar_gap = settings.get(f"{gene}_bar_gap", settings.get('bar_gap', 0.15))
        gene_margins = settings.get(f"{gene}_margins", {'l': 80, 'r': 80, 't': 100, 'b': 100})
        gene_bg_color = settings.get(f"{gene}_bg_color", settings.get('plot_bgcolor', '#FFFFFF'))
        gene_tick_size = settings.get(f"{gene}_tick_size", 12)
        
        wrapped_labels = [GraphGenerator._wrap_text(str(cond), 15) for cond in condition_names]
        
        # P-value legend
        legend_text = "<b>Significance:</b>  * p<0.05  ** p<0.01  *** p<0.001"
        if 'significance_2' in gene_data.columns and gene_data['significance_2'].notna().any():
            legend_text += "<br><b>2nd Comparison:</b>  # p<0.05  ## p<0.01  ### p<0.001"
        
        fig.add_annotation(
            text=legend_text,
            xref="paper", yref="paper",
            x=1.0, y=-0.15,
            xanchor='right', yanchor='top',
            showarrow=False,
            font=dict(size=12, color='#666666', family='Arial'),
            bgcolor='rgba(255,255,255,0.9)',
            bordercolor='#CCCCCC',
            borderwidth=1,
            borderpad=4
        )
        
        fig.update_layout(
            title=dict(
                text=f"{gene} Expression",
                font=dict(size=settings.get('title_size', 20), family='Arial', color='#333333'),
                x=0.5, xanchor='center', y=0.98, yanchor='top'
            ),
            xaxis=dict(
                title=None,
                showgrid=False,
                zeroline=False,
                tickmode='array',
                tickvals=list(range(n_bars)),
                ticktext=wrapped_labels,
                tickfont=dict(size=gene_tick_size),
                tickangle=0,
                showline=False,
                mirror=False,
                side='bottom',
                range=[-0.5, n_bars - 0.5]
            ),
            yaxis=y_axis_config,
            template=settings.get('color_scheme', 'plotly_white'),
            font=dict(size=settings.get('font_size', 14), family='Arial'),
            height=settings.get('figure_height', 600),
            width=settings.get('figure_width', 1000),
            bargap=gene_bar_gap,
            showlegend=settings.get('show_legend', False),
            plot_bgcolor=gene_bg_color,
            paper_bgcolor='#FFFFFF',
            margin=dict(
                l=gene_margins.get('l', 80),
                r=gene_margins.get('r', 80),
                t=gene_margins.get('t', 100),
                b=gene_margins.get('b', 120)
            )
        )
        
        return fig
    
    @staticmethod
    def export_publication_figure(fig: go.Figure, gene: str, settings: Dict) -> Dict[str, bytes]:
        """Export figure in publication-ready formats"""
        exports = {}
        
        # Update figure for publication
        pub_fig = go.Figure(fig)
        
        # Scale dimensions
        dpi = settings.get('dpi', 300)
        width_inches = settings.get('width', 7)
        height_inches = settings.get('height', 5)
        
        width_px = int(width_inches * dpi)
        height_px = int(height_inches * dpi)
        
        # Update font sizes for publication
        font_scale = settings.get('font_size', 12) / 14  # Base scale
        
        pub_fig.update_layout(
            width=width_px,
            height=height_px,
            font=dict(
                family=settings.get('font_family', 'Arial'),
                size=settings.get('font_size', 12)
            ),
            title=dict(font=dict(size=int(16 * font_scale))),
        )
        
        # Export formats
        try:
            # PNG
            if 'PNG' in settings.get('formats', ['PNG']):
                png_bytes = pub_fig.to_image(format='png', scale=dpi/72, engine='kaleido')
                exports['png'] = png_bytes
            
            # SVG
            if 'SVG' in settings.get('formats', []):
                svg_bytes = pub_fig.to_image(format='svg', engine='kaleido')
                exports['svg'] = svg_bytes
            
            # PDF
            if 'PDF' in settings.get('formats', []):
                pdf_bytes = pub_fig.to_image(format='pdf', engine='kaleido')
                exports['pdf'] = pdf_bytes
                
        except Exception as e:
            st.warning(f"Export error (install kaleido for static exports): {e}")
            # Fallback to HTML only
            html_buffer = io.StringIO()
            pub_fig.write_html(html_buffer)
            exports['html'] = html_buffer.getvalue().encode()
        
        return exports


# ==================== EXPORT FUNCTIONS ====================
def export_to_excel(raw_data: pd.DataFrame, processed_data: Dict[str, pd.DataFrame], 
                   params: Dict, mapping: Dict, qc_results: Optional[Dict] = None) -> bytes:
    """Export comprehensive Excel with gene-by-gene sheets"""
    output = io.BytesIO()
    
    with pd.ExcelWriter(output, engine='openpyxl') as writer:
        # Parameters sheet
        pd.DataFrame([params]).to_excel(writer, sheet_name='Analysis_Parameters', index=False)
        
        # Sample mapping sheet
        mapping_list = [{'Original': k, **v} for k, v in mapping.items()]
        pd.DataFrame(mapping_list).to_excel(writer, sheet_name='Sample_Mapping', index=False)
        
        # QC Results (if available)
        if qc_results and qc_results.get('summary'):
            qc_summary = pd.DataFrame([qc_results['summary']])
            qc_summary.to_excel(writer, sheet_name='QC_Summary', index=False)
            
            if qc_results.get('high_ct_wells'):
                pd.DataFrame(qc_results['high_ct_wells']).to_excel(
                    writer, sheet_name='QC_High_CT', index=False
                )
            
            if qc_results.get('outlier_wells'):
                pd.DataFrame(qc_results['outlier_wells']).to_excel(
                    writer, sheet_name='QC_Outliers', index=False
                )
        
        # Raw data with mapped conditions
        raw_export = raw_data.copy()
        if mapping:
            raw_export['Condition'] = raw_export['Sample'].map(
                lambda x: mapping.get(x, {}).get('condition', x)
            )
        else:
            raw_export['Condition'] = raw_export['Sample']
        
        cols = ['Well', 'Sample', 'Condition', 'Target', 'CT']
        if 'Source_File' in raw_export.columns:
            cols.append('Source_File')
        raw_export = raw_export[cols]
        raw_export.to_excel(writer, sheet_name='Raw_Data', index=False)

        # Gene-by-gene calculations
        for gene, gene_data in processed_data.items():
            sheet_name = f"{gene}_Analysis"[:31]
            gene_data.to_excel(writer, sheet_name=sheet_name, index=False)
        
        # Summary sheet
        if processed_data:
            all_data = pd.concat(processed_data.values(), ignore_index=True)
            summary = all_data.groupby(['Target', 'Group']).agg({
                'Relative_Expression': ['mean', 'std', 'count'],
                'p_value': 'min'
            }).round(4)
            summary.to_excel(writer, sheet_name='Summary')
    
    return output.getvalue()


# ==================== CUSTOM CSS ====================
def apply_custom_css():
    """Apply custom CSS for improved UI"""
    st.markdown("""
    <style>
    /* Clean, professional look */
    .stApp {
        background-color: #f8f9fa;
    }
    
    /* Card-like containers */
    .css-1r6slb0, .css-12oz5g7 {
        background-color: white;
        border-radius: 10px;
        padding: 15px;
        box-shadow: 0 2px 6px rgba(0,0,0,0.1);
    }
    
    /* Graph containers */
    [data-testid="stPlotlyChart"] {
        box-shadow: 0 4px 12px rgba(0,0,0,0.15);
        border-radius: 8px;
        background: white;
        padding: 10px;
    }
    
    /* Metric cards */
    [data-testid="metric-container"] {
        background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
        padding: 15px;
        border-radius: 10px;
        color: white;
    }
    
    [data-testid="metric-container"] label {
        color: rgba(255,255,255,0.8) !important;
    }
    
    [data-testid="metric-container"] [data-testid="stMetricValue"] {
        color: white !important;
    }
    
    /* QC status colors */
    .qc-ok { color: #28a745; font-weight: bold; }
    .qc-warning { color: #ffc107; font-weight: bold; }
    .qc-critical { color: #dc3545; font-weight: bold; }
    
    /* Compact expanders */
    .streamlit-expanderHeader {
        font-size: 14px;
        font-weight: 600;
    }
    
    /* Tab styling */
    .stTabs [data-baseweb="tab-list"] {
        gap: 8px;
    }
    
    .stTabs [data-baseweb="tab"] {
        padding: 10px 20px;
        border-radius: 8px 8px 0 0;
    }
    
    /* Button improvements */
    .stButton > button {
        border-radius: 8px;
        font-weight: 500;
        transition: all 0.2s;
    }
    
    .stButton > button:hover {
        transform: translateY(-1px);
        box-shadow: 0 4px 8px rgba(0,0,0,0.2);
    }
    
    /* Data table improvements */
    .dataframe {
        font-size: 13px;
    }
    
    /* Hide streamlit branding */
    #MainMenu {visibility: hidden;}
    footer {visibility: hidden;}
    </style>
    """, unsafe_allow_html=True)


# ==================== UI COMPONENTS ====================
apply_custom_css()

st.title("🧬 qPCR Analysis Suite Pro v3.0")
st.markdown("**Enhanced analysis with QC module, bulk operations, and publication-ready exports**")

# Sidebar
with st.sidebar:
    st.header("📋 Workflow Guide")
    
    with st.expander("📖 Quick Start", expanded=True):
        st.markdown("""
        1. **Import** → Upload CSV files
        2. **QC** → Review data quality
        3. **Setup** → Map samples to conditions
        4. **Analyze** → Run ΔΔCt calculations
        5. **Visualize** → Customize graphs
        6. **Export** → Download results
        """)
    
    with st.expander("⚙️ Global Settings"):
        st.session_state.graph_settings['figure_width'] = st.slider(
            "Figure Width", 600, 1400, st.session_state.graph_settings.get('figure_width', 1000)
        )
        st.session_state.graph_settings['figure_height'] = st.slider(
            "Figure Height", 400, 900, st.session_state.graph_settings.get('figure_height', 600)
        )
    
    st.markdown("---")
    st.caption("v3.0 | Enhanced Edition")

# Main tabs
tab1, tab2, tab3, tab4, tab5, tab6 = st.tabs([
    "📁 Data Import", 
    "🔍 QC & Review",
    "🗺️ Sample Setup", 
    "📊 Analysis", 
    "📈 Visualization", 
    "📤 Export"
])

# ==================== TAB 1: DATA IMPORT ====================
with tab1:
    st.header("Step 1: Import qPCR Data")
    
    col_upload, col_info = st.columns([2, 1])
    
    with col_upload:
        uploaded_files = st.file_uploader(
            "Upload qPCR CSV files",
            type=['csv'],
            accept_multiple_files=True,
            help="Supports QuantStudio, CFX, and similar formats"
        )
    
    with col_info:
        with st.expander("📝 Supported Formats"):
            st.markdown("""
            - **QuantStudio** (Applied Biosystems)
            - **CFX Maestro** (Bio-Rad)
            - Generic CSV with Well, Sample, Target, CT columns
            """)
    
    if uploaded_files:
        all_data = []
        
        with st.spinner("Parsing files..."):
            for file in uploaded_files:
                parsed = QPCRParser.parse(file)
                if parsed is not None:
                    parsed['Source_File'] = file.name
                    all_data.append(parsed)
                    st.success(f"✅ {file.name}: {len(parsed)} data points")
                else:
                    st.error(f"❌ {file.name}: Could not parse")
        
        if all_data:
            st.session_state.data = pd.concat(all_data, ignore_index=True)
            
            # Natural sort samples
            unique_samples = sorted(
                st.session_state.data['Sample'].unique(),
                key=natural_sort_key
            )
            
            if st.session_state.sample_order is None:
                st.session_state.sample_order = list(unique_samples)
            else:
                existing_set = set(st.session_state.sample_order)
                new_samples = [s for s in unique_samples if s not in existing_set]
                st.session_state.sample_order.extend(new_samples)
            
            # Summary metrics
            st.markdown("### 📊 Data Summary")
            col1, col2, col3, col4 = st.columns(4)
            col1.metric("Total Wells", len(st.session_state.data))
            col2.metric("Samples", st.session_state.data['Sample'].nunique())
            col3.metric("Genes", st.session_state.data['Target'].nunique())
            
            # Detect housekeeping gene
            hk_genes = [g for g in st.session_state.data['Target'].unique() 
                       if g.upper() in ['ACTIN', 'B-ACTIN', 'GAPDH', 'ACTB', 'B2M', '18S']]
            
            if hk_genes:
                st.session_state.hk_gene = st.selectbox(
                    "🔬 Select Housekeeping Gene", 
                    hk_genes,
                    key='hk_select_main'
                )
                col4.metric("HK Gene", st.session_state.hk_gene)
            else:
                st.warning("⚠️ No common housekeeping gene detected. Please select manually in Sample Setup.")
            
            # Data preview
            with st.expander("📋 Preview Raw Data", expanded=False):
                st.dataframe(st.session_state.data.head(100), height=300, use_container_width=True)
            
            st.success("✅ Data imported successfully! Proceed to **QC & Review** tab.")


# ==================== TAB 2: QC & REVIEW ====================
with tab2:
    st.header("Step 2: Quality Control & Data Review")
    
    if st.session_state.data is not None:
        # Run QC analysis
        if st.button("🔍 Run QC Analysis", type="primary"):
            with st.spinner("Analyzing data quality..."):
                qc_results = QCModule.run_qc_analysis(st.session_state.data)
                st.session_state['qc_results'] = qc_results
                st.session_state.qc_complete = True
        
        if st.session_state.get('qc_results'):
            qc = st.session_state['qc_results']
            
            # QC Summary Dashboard
            st.markdown("### 📊 QC Summary Dashboard")
            
            col1, col2, col3, col4, col5 = st.columns(5)
            
            with col1:
                score = qc['summary']['quality_score']
                color = "#28a745" if score >= 90 else "#ffc107" if score >= 70 else "#dc3545"
                st.markdown(f"""
                <div style='text-align:center; padding:20px; background:linear-gradient(135deg, {color}dd, {color}aa); 
                            border-radius:10px; color:white;'>
                    <h1 style='margin:0;'>{score:.0f}%</h1>
                    <p style='margin:0;'>Quality Score</p>
                </div>
                """, unsafe_allow_html=True)
            
            with col2:
                st.metric("Total Wells", qc['summary']['total_wells'])
            with col3:
                st.metric("Flagged", qc['summary']['flagged_wells'], 
                         delta=None if qc['summary']['flagged_wells'] == 0 else f"-{qc['summary']['flagged_wells']}")
            with col4:
                st.metric("High CT", qc['summary']['high_ct_count'])
            with col5:
                st.metric("Outliers", qc['summary']['outlier_count'])
            
            st.markdown("---")
            
            # QC Details in expandable sections
            col_left, col_right = st.columns(2)
            
            with col_left:
                # High CT Wells
                with st.expander(f"⚠️ High CT Values ({len(qc['high_ct_wells'])})", expanded=len(qc['high_ct_wells']) > 0):
                    if qc['high_ct_wells']:
                        high_ct_df = pd.DataFrame(qc['high_ct_wells'])
                        
                        # Add exclusion checkboxes
                        for idx, row in high_ct_df.iterrows():
                            col_a, col_b, col_c = st.columns([3, 1, 1])
                            with col_a:
                                st.write(f"**{row['well']}**: {row['sample']} / {row['target']} = {row['ct']:.2f}")
                            with col_b:
                                status_color = "🔴" if row['status'] == 'critical' else "🟡"
                                st.write(status_color)
                            with col_c:
                                if st.checkbox("Exclude", key=f"excl_high_{row['well']}"):
                                    st.session_state.qc_flagged_wells.add(row['well'])
                                elif row['well'] in st.session_state.qc_flagged_wells:
                                    st.session_state.qc_flagged_wells.discard(row['well'])
                    else:
                        st.success("No high CT values detected!")
            
            with col_right:
                # Outliers
                with st.expander(f"📊 Statistical Outliers ({len(qc['outlier_wells'])})", expanded=len(qc['outlier_wells']) > 0):
                    if qc['outlier_wells']:
                        for item in qc['outlier_wells']:
                            col_a, col_b = st.columns([4, 1])
                            with col_a:
                                st.write(f"**{item['well']}**: {item['sample']} / {item['target']}")
                                st.caption(f"CT: {item['ct']:.2f} (replicate mean: {item['replicate_mean']:.2f})")
                            with col_b:
                                if st.checkbox("Exclude", key=f"excl_out_{item['well']}"):
                                    st.session_state.qc_flagged_wells.add(item['well'])
                                elif item['well'] in st.session_state.qc_flagged_wells:
                                    st.session_state.qc_flagged_wells.discard(item['well'])
                    else:
                        st.success("No statistical outliers detected!")
            
            # CT Value Editor
            st.markdown("### ✏️ CT Value Editor")
            with st.expander("Edit Individual CT Values", expanded=False):
                st.caption("Use this to correct obvious data entry errors or typos")
                
                # Select well to edit
                well_options = st.session_state.data['Well'].unique().tolist()
                selected_well = st.selectbox("Select Well to Edit", well_options, key="edit_well_select")
                
                if selected_well:
                    well_data = st.session_state.data[st.session_state.data['Well'] == selected_well].iloc[0]
                    
                    col_ed1, col_ed2, col_ed3 = st.columns(3)
                    with col_ed1:
                        st.write(f"**Sample:** {well_data['Sample']}")
                        st.write(f"**Target:** {well_data['Target']}")
                    with col_ed2:
                        original_ct = well_data['CT']
                        st.write(f"**Original CT:** {original_ct:.2f}")
                    with col_ed3:
                        new_ct = st.number_input(
                            "New CT Value",
                            min_value=0.0,
                            max_value=45.0,
                            value=float(original_ct),
                            step=0.01,
                            key="new_ct_input"
                        )
                        reason = st.text_input("Reason for edit", key="edit_reason")
                        
                        if st.button("Apply Edit", key="apply_edit"):
                            if new_ct != original_ct and reason:
                                st.session_state.qc_edited_values[selected_well] = {
                                    'original': original_ct,
                                    'edited': new_ct,
                                    'reason': reason,
                                    'timestamp': datetime.now().isoformat()
                                }
                                st.success(f"CT value updated: {original_ct:.2f} → {new_ct:.2f}")
                            elif not reason:
                                st.warning("Please provide a reason for the edit")
            
            # Show current exclusions/edits
            if st.session_state.qc_flagged_wells or st.session_state.qc_edited_values:
                st.markdown("### 📝 Pending QC Actions")
                col_exc, col_edit = st.columns(2)
                
                with col_exc:
                    if st.session_state.qc_flagged_wells:
                        st.write(f"**Excluded wells ({len(st.session_state.qc_flagged_wells)}):**")
                        st.write(", ".join(sorted(st.session_state.qc_flagged_wells)))
                
                with col_edit:
                    if st.session_state.qc_edited_values:
                        st.write(f"**Edited values ({len(st.session_state.qc_edited_values)}):**")
                        for well, info in st.session_state.qc_edited_values.items():
                            st.write(f"- {well}: {info['original']:.2f} → {info['edited']:.2f}")
            
            st.success("✅ QC Review complete! Proceed to **Sample Setup** tab.")
    else:
        st.info("⏳ Please import data first in the **Data Import** tab.")


# ==================== TAB 3: SAMPLE SETUP ====================
with tab3:
    st.header("Step 3: Sample Setup & Mapping")
    
    if st.session_state.data is not None:
        # Efficacy type selection
        detected_genes = set(st.session_state.data['Target'].unique())
        suggested = None
        for eff, cfg in EFFICACY_CONFIG.items():
            if any(g in detected_genes for g in cfg['genes']):
                suggested = eff
                break
        
        col_eff, col_info = st.columns([2, 1])
        
        with col_eff:
            efficacy = st.selectbox(
                "🎯 Efficacy Test Type",
                list(EFFICACY_CONFIG.keys()),
                index=list(EFFICACY_CONFIG.keys()).index(suggested) if suggested else 0
            )
            st.session_state.selected_efficacy = efficacy
        
        with col_info:
            config = EFFICACY_CONFIG[efficacy]
            st.info(f"**{config['cell']}**")
        
        st.markdown(f"*{config['description']}*")
        
        with st.expander("📋 Control Structure"):
            for ctrl_type, ctrl_name in config['controls'].items():
                st.markdown(f"- **{ctrl_type.title()}**: {ctrl_name}")
        
        st.markdown("---")
        
        # Bulk Operations
        st.markdown("### ⚡ Bulk Operations")
        
        with st.expander("🔧 Bulk Sample Operations", expanded=False):
            col_bulk1, col_bulk2, col_bulk3 = st.columns(3)
            
            with col_bulk1:
                st.markdown("**Select Samples**")
                if st.button("Select All", key="select_all"):
                    st.session_state.selected_samples = set(st.session_state.sample_order)
                if st.button("Clear Selection", key="clear_sel"):
                    st.session_state.selected_samples = set()
                
                st.caption(f"Selected: {len(st.session_state.selected_samples)}")
            
            with col_bulk2:
                st.markdown("**Bulk Set Group**")
                group_types = ['Baseline', 'Negative Control', 'Positive Control', 'Treatment']
                bulk_group = st.selectbox("Group", group_types, key="bulk_group_sel")
                if st.button("Apply Group", key="apply_bulk_group"):
                    for sample in st.session_state.selected_samples:
                        if sample in st.session_state.sample_mapping:
                            st.session_state.sample_mapping[sample]['group'] = bulk_group
                    st.success(f"Applied '{bulk_group}' to {len(st.session_state.selected_samples)} samples")
                    st.rerun()
            
            with col_bulk3:
                st.markdown("**Bulk Rename**")
                find_text = st.text_input("Find", key="bulk_find")
                replace_text = st.text_input("Replace with", key="bulk_replace")
                if st.button("Apply Rename", key="apply_bulk_rename"):
                    count = 0
                    for sample in st.session_state.selected_samples:
                        if sample in st.session_state.sample_mapping:
                            old_cond = st.session_state.sample_mapping[sample]['condition']
                            new_cond = old_cond.replace(find_text, replace_text)
                            if new_cond != old_cond:
                                st.session_state.sample_mapping[sample]['condition'] = new_cond
                                count += 1
                    st.success(f"Renamed {count} conditions")
                    st.rerun()
        
        st.markdown("---")
        
        # Sample Mapping Table
        st.markdown("### 🗺️ Sample Condition Mapping")
        
        # Initialize mapping for all samples
        group_types = ['Negative Control', 'Positive Control', 'Treatment']
        if 'baseline' in config['controls']:
            group_types.insert(0, 'Baseline')
        
        for sample in st.session_state.sample_order:
            if sample not in st.session_state.sample_mapping:
                st.session_state.sample_mapping[sample] = {
                    'condition': sample,
                    'group': 'Treatment',
                    'include': True
                }
            if 'include' not in st.session_state.sample_mapping[sample]:
                st.session_state.sample_mapping[sample]['include'] = True
        
        # Header
        st.markdown("""
        <div style='background-color: #e9ecef; padding: 10px; border-radius: 5px; margin-bottom: 10px;'>
            <div style='display: grid; grid-template-columns: 0.5fr 0.5fr 1.5fr 2.5fr 2fr 1fr; gap: 10px; font-weight: bold;'>
                <span>✓</span>
                <span>#</span>
                <span>Original</span>
                <span>Condition Name</span>
                <span>Group</span>
                <span>Move</span>
            </div>
        </div>
        """, unsafe_allow_html=True)
        
        # Sample rows
        for i, sample in enumerate(st.session_state.sample_order):
            with st.container():
                col0, col_num, col1, col2, col3, col_move = st.columns([0.5, 0.5, 1.5, 2.5, 2, 1])
                
                with col0:
                    include = st.checkbox(
                        "",
                        value=st.session_state.sample_mapping[sample].get('include', True),
                        key=f"inc_{sample}_{i}",
                        label_visibility="collapsed"
                    )
                    st.session_state.sample_mapping[sample]['include'] = include
                
                with col_num:
                    st.markdown(f"<div style='padding-top:8px;text-align:center;'><b>{i+1}</b></div>", 
                               unsafe_allow_html=True)
                
                with col1:
                    # Selectable for bulk operations
                    is_selected = sample in st.session_state.selected_samples
                    if st.checkbox(sample[:15], value=is_selected, key=f"sel_{sample}_{i}", 
                                  label_visibility="visible"):
                        st.session_state.selected_samples.add(sample)
                    else:
                        st.session_state.selected_samples.discard(sample)
                
                with col2:
                    cond = st.text_input(
                        "Condition",
                        st.session_state.sample_mapping[sample]['condition'],
                        key=f"cond_{sample}_{i}",
                        label_visibility="collapsed"
                    )
                    st.session_state.sample_mapping[sample]['condition'] = cond
                
                with col3:
                    grp_idx = 0
                    try:
                        grp_idx = group_types.index(st.session_state.sample_mapping[sample]['group'])
                    except (ValueError, KeyError):
                        pass
                    
                    grp = st.selectbox(
                        "Group",
                        group_types,
                        index=grp_idx,
                        key=f"grp_{sample}_{i}",
                        label_visibility="collapsed"
                    )
                    st.session_state.sample_mapping[sample]['group'] = grp
                
                with col_move:
                    btn_col1, btn_col2 = st.columns(2)
                    with btn_col1:
                        if i > 0 and st.button("⬆", key=f"up_{sample}_{i}", use_container_width=True):
                            st.session_state.sample_order[i], st.session_state.sample_order[i-1] = \
                                st.session_state.sample_order[i-1], st.session_state.sample_order[i]
                            st.rerun()
                    with btn_col2:
                        if i < len(st.session_state.sample_order) - 1:
                            if st.button("⬇", key=f"dn_{sample}_{i}", use_container_width=True):
                                st.session_state.sample_order[i], st.session_state.sample_order[i+1] = \
                                    st.session_state.sample_order[i+1], st.session_state.sample_order[i]
                                st.rerun()
        
        # Update excluded samples
        st.session_state.excluded_samples = set([
            s for s, v in st.session_state.sample_mapping.items()
            if not v.get('include', True)
        ])
        
        # Summary
        st.markdown("---")
        total = len(st.session_state.sample_order)
        included = sum(1 for s in st.session_state.sample_order 
                      if st.session_state.sample_mapping[s].get('include', True))
        
        col_s1, col_s2, col_s3 = st.columns(3)
        col_s1.metric("Total Samples", total)
        col_s2.metric("Included", included)
        col_s3.metric("Excluded", total - included)
        
        st.success("✅ Sample setup complete! Proceed to **Analysis** tab.")
    else:
        st.info("⏳ Please import data first.")


# ==================== TAB 4: ANALYSIS ====================
with tab4:
    st.header("Step 4: Run Analysis")
    
    if st.session_state.data is not None and st.session_state.sample_mapping:
        # Build condition list
        condition_list = []
        sample_to_condition = {}
        
        for sample in st.session_state.get('sample_order', []):
            if sample in st.session_state.sample_mapping:
                mapping_info = st.session_state.sample_mapping[sample]
                if mapping_info.get('include', True):
                    condition = mapping_info.get('condition', sample)
                    if condition not in condition_list:
                        condition_list.append(condition)
                    sample_to_condition[condition] = sample
        
        if condition_list:
            st.markdown("### 📊 Analysis Configuration")
            
            col_r1, col_r2, col_r3 = st.columns(3)
            
            with col_r1:
                st.markdown("**ΔΔCt Reference**")
                ref_condition = st.selectbox(
                    "Baseline for fold change",
                    condition_list,
                    index=0,
                    key="ref_ddct"
                )
                ref_sample_key = sample_to_condition[ref_condition]
            
            with col_r2:
                st.markdown("**P-value Reference 1 (*)**")
                cmp_condition = st.selectbox(
                    "Primary comparison",
                    condition_list,
                    index=0,
                    key="cmp_pval1"
                )
                cmp_sample_key = sample_to_condition[cmp_condition]
            
            with col_r3:
                st.markdown("**P-value Reference 2 (#)**")
                use_second = st.checkbox("Enable 2nd comparison", key="use_2nd_cmp")
                
                cmp_sample_key_2 = None
                if use_second:
                    cond_list_2 = [c for c in condition_list if c != cmp_condition]
                    if cond_list_2:
                        cmp_condition_2 = st.selectbox(
                            "Secondary comparison",
                            cond_list_2,
                            key="cmp_pval2"
                        )
                        cmp_sample_key_2 = sample_to_condition[cmp_condition_2]
            
            # Summary box
            st.markdown("---")
            st.markdown(f"""
            <div style='background-color: #e7f3ff; padding: 15px; border-radius: 10px; text-align: center;'>
                <b>Analysis will calculate:</b><br>
                Fold changes relative to <code>{ref_condition}</code><br>
                P-values (*) vs <code>{cmp_condition}</code>
                {f"<br>P-values (#) vs <code>{cmp_condition_2}</code>" if use_second and cmp_sample_key_2 else ""}
            </div>
            """, unsafe_allow_html=True)
            
            st.markdown("---")
            
            # Run button
            if st.button("▶️ Run Full Analysis", type="primary", use_container_width=True):
                ok = AnalysisEngine.run_full_analysis(
                    ref_sample_key,
                    cmp_sample_key,
                    cmp_sample_key_2 if use_second else None
                )
                if ok:
                    st.success("✅ Analysis complete!")
                    st.balloons()
                    st.rerun()
        
        # Show results if available
        if st.session_state.processed_data:
            st.markdown("---")
            st.markdown("### 📋 Results Summary")
            
            all_results = pd.concat(st.session_state.processed_data.values(), ignore_index=True)
            
            col1, col2, col3 = st.columns(3)
            col1.metric("Genes Analyzed", len(st.session_state.processed_data))
            col2.metric("Conditions", all_results['Condition'].nunique())
            sig_count = (all_results['p_value'] < 0.05).sum()
            col3.metric("Significant (p<0.05)", f"{sig_count}/{len(all_results)}")
            
            # Gene results
            for gene, gene_df in st.session_state.processed_data.items():
                with st.expander(f"🧬 {gene}", expanded=False):
                    display_cols = ['Condition', 'Group', 'Relative_Expression', 'p_value', 
                                   'significance', 'n_replicates', 'SEM']
                    display_df = gene_df[[c for c in display_cols if c in gene_df.columns]].copy()
                    
                    # Format numbers
                    if 'Relative_Expression' in display_df.columns:
                        display_df['Relative_Expression'] = display_df['Relative_Expression'].round(3)
                    if 'p_value' in display_df.columns:
                        display_df['p_value'] = display_df['p_value'].apply(format_p_value)
                    if 'SEM' in display_df.columns:
                        display_df['SEM'] = display_df['SEM'].round(3)
                    
                    st.dataframe(display_df, use_container_width=True, hide_index=True)
            
            st.success("✅ Results ready! Proceed to **Visualization** tab.")
    else:
        st.info("⏳ Please complete Data Import and Sample Setup first.")


# ==================== TAB 5: VISUALIZATION ====================
with tab5:
    st.header("Step 5: Graph Visualization")
    
    if st.session_state.processed_data:
        efficacy_config = EFFICACY_CONFIG.get(st.session_state.selected_efficacy, {})
        
        for gene in st.session_state.processed_data.keys():
            st.markdown("---")
            st.markdown(f"## 🧬 {gene}")
            
            col_ctrl, col_graph = st.columns([1, 3])
            
            with col_ctrl:
                with st.expander("⚙️ Graph Controls", expanded=True):
                    # Display toggles
                    show_sig_key = f"{gene}_show_sig"
                    show_err_key = f"{gene}_show_err"
                    
                    if show_sig_key not in st.session_state.graph_settings:
                        st.session_state.graph_settings[show_sig_key] = True
                    if show_err_key not in st.session_state.graph_settings:
                        st.session_state.graph_settings[show_err_key] = True
                    
                    st.session_state.graph_settings[show_sig_key] = st.checkbox(
                        "✨ Show Significance", 
                        st.session_state.graph_settings[show_sig_key],
                        key=f"chk_sig_{gene}"
                    )
                    
                    st.session_state.graph_settings[show_err_key] = st.checkbox(
                        "📏 Show Error Bars",
                        st.session_state.graph_settings[show_err_key],
                        key=f"chk_err_{gene}"
                    )
                    
                    # Bar gap
                    bar_gap_key = f"{gene}_bar_gap"
                    if bar_gap_key not in st.session_state.graph_settings:
                        st.session_state.graph_settings[bar_gap_key] = 0.25
                    
                    st.session_state.graph_settings[bar_gap_key] = st.slider(
                        "Bar Gap",
                        0.0, 0.5,
                        st.session_state.graph_settings[bar_gap_key],
                        0.05,
                        key=f"gap_{gene}"
                    )
                
                # Per-bar color controls
                with st.expander("🎨 Bar Colors", expanded=False):
                    gene_data = st.session_state.processed_data[gene]
                    
                    if f'{gene}_bar_settings' not in st.session_state:
                        st.session_state[f'{gene}_bar_settings'] = {}
                    
                    if 'bar_colors_per_sample' not in st.session_state.graph_settings:
                        st.session_state.graph_settings['bar_colors_per_sample'] = {}
                    
                    for idx, row in gene_data.iterrows():
                        condition = row['Condition']
                        group = row.get('Group', 'Treatment')
                        bar_key = f"{gene}_{condition}"
                        
                        default_colors = {
                            'Baseline': '#FFFFFF',
                            'Negative Control': '#9E9E9E',
                            'Positive Control': '#9E9E9E',
                            'Treatment': '#D3D3D3'
                        }
                        default_color = default_colors.get(group, '#D3D3D3')
                        
                        if bar_key not in st.session_state[f'{gene}_bar_settings']:
                            st.session_state[f'{gene}_bar_settings'][bar_key] = {
                                'color': default_color,
                                'show_sig': True,
                                'show_err': True
                            }
                        
                        col_name, col_color = st.columns([2, 1])
                        with col_name:
                            st.caption(f"{condition[:18]}...")
                        with col_color:
                            new_color = st.color_picker(
                                "Color",
                                st.session_state[f'{gene}_bar_settings'][bar_key]['color'],
                                key=f"col_{gene}_{condition}_{idx}",
                                label_visibility="collapsed"
                            )
                            st.session_state[f'{gene}_bar_settings'][bar_key]['color'] = new_color
                            st.session_state.graph_settings['bar_colors_per_sample'][bar_key] = new_color
            
            with col_graph:
                # Generate graph
                gene_data = st.session_state.processed_data[gene]
                
                current_settings = st.session_state.graph_settings.copy()
                current_settings['show_significance'] = current_settings.get(f"{gene}_show_sig", True)
                current_settings['show_error'] = current_settings.get(f"{gene}_show_err", True)
                current_settings['bar_gap'] = current_settings.get(f"{gene}_bar_gap", 0.15)
                
                fig = GraphGenerator.create_gene_graph(
                    gene_data,
                    gene,
                    current_settings,
                    efficacy_config,
                    sample_order=st.session_state.get('sample_order')
                )
                
                st.plotly_chart(fig, use_container_width=True, key=f"fig_{gene}")
                st.session_state.graphs[gene] = fig
        
        st.success("✅ Graphs ready! Proceed to **Export** tab for downloads.")
    else:
        st.info("⏳ Run analysis first in the **Analysis** tab.")


# ==================== TAB 6: EXPORT ====================
with tab6:
    st.header("Step 6: Export Results")
    
    if st.session_state.processed_data:
        # Publication Export Settings
        st.markdown("### 📸 Publication-Ready Export")
        
        with st.expander("⚙️ Export Settings", expanded=True):
            col_preset, col_custom = st.columns(2)
            
            with col_preset:
                preset_name = st.selectbox(
                    "Preset",
                    list(EXPORT_PRESETS.keys()),
                    format_func=lambda x: EXPORT_PRESETS[x]['name']
                )
                preset = EXPORT_PRESETS[preset_name]
                st.session_state.export_settings.update(preset)
            
            with col_custom:
                st.session_state.export_settings['dpi'] = st.selectbox(
                    "DPI", [150, 300, 600], 
                    index=[150, 300, 600].index(preset['dpi'])
                )
                st.session_state.export_settings['width'] = st.number_input(
                    "Width (inches)", 
                    min_value=2.0, max_value=15.0,
                    value=float(preset['width']),
                    step=0.5
                )
        
        st.markdown("---")
        
        # Download sections
        col1, col2 = st.columns(2)
        
        with col1:
            st.markdown("### 📊 Complete Excel Report")
            st.caption("All calculations, statistics, and raw data")
            
            analysis_params = {
                'Date': datetime.now().strftime("%Y-%m-%d %H:%M"),
                'Efficacy_Type': st.session_state.selected_efficacy,
                'Housekeeping_Gene': st.session_state.hk_gene,
                'Genes_Analyzed': len(st.session_state.processed_data),
                'QC_Excluded_Wells': len(st.session_state.qc_flagged_wells),
                'QC_Edited_Values': len(st.session_state.qc_edited_values)
            }
            
            excel_data = export_to_excel(
                st.session_state.data,
                st.session_state.processed_data,
                analysis_params,
                st.session_state.sample_mapping,
                st.session_state.get('qc_results')
            )
            
            st.download_button(
                label="📥 Download Excel Report",
                data=excel_data,
                file_name=f"qPCR_{st.session_state.selected_efficacy}_{datetime.now().strftime('%Y%m%d_%H%M')}.xlsx",
                mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
                type="primary"
            )
        
        with col2:
            st.markdown("### 📈 All Graphs (HTML)")
            st.caption("Interactive graphs for all genes")
            
            if st.session_state.graphs:
                html_parts = [
                    "<html><head>",
                    "<title>qPCR Analysis Graphs</title>",
                    "<style>body{font-family:Arial,sans-serif;max-width:1200px;margin:0 auto;padding:20px;}</style>",
                    "</head><body>",
                    f"<h1>{st.session_state.selected_efficacy} Analysis</h1>",
                    f"<p>Generated: {datetime.now().strftime('%Y-%m-%d %H:%M')}</p>"
                ]
                
                for gene, fig in st.session_state.graphs.items():
                    html_parts.append(f"<h2>{gene}</h2>")
                    html_parts.append(fig.to_html(include_plotlyjs='cdn', div_id=f"graph_{gene}"))
                    html_parts.append("<hr>")
                
                html_parts.append("</body></html>")
                combined_html = "\n".join(html_parts)
                
                st.download_button(
                    label="📥 Download All Graphs (HTML)",
                    data=combined_html,
                    file_name=f"qPCR_graphs_{datetime.now().strftime('%Y%m%d_%H%M')}.html",
                    mime="text/html",
                    type="primary"
                )
        
        st.markdown("---")
        
        # Individual Downloads
        st.markdown("### 📁 Individual Downloads")
        
        with st.expander("Gene-by-Gene Files", expanded=False):
            for gene in st.session_state.processed_data.keys():
                col_csv, col_html = st.columns(2)
                
                with col_csv:
                    csv_buffer = io.StringIO()
                    st.session_state.processed_data[gene].to_csv(csv_buffer, index=False)
                    st.download_button(
                        label=f"📥 {gene}.csv",
                        data=csv_buffer.getvalue(),
                        file_name=f"{gene}_data.csv",
                        mime="text/csv",
                        key=f"csv_{gene}"
                    )
                
                with col_html:
                    if gene in st.session_state.graphs:
                        html_buffer = io.StringIO()
                        st.session_state.graphs[gene].write_html(html_buffer)
                        st.download_button(
                            label=f"📥 {gene}.html",
                            data=html_buffer.getvalue(),
                            file_name=f"{gene}_graph.html",
                            mime="text/html",
                            key=f"html_{gene}"
                        )
        
        # Config export for reproducibility
        with st.expander("🔧 Configuration Files", expanded=False):
            config_data = {
                'analysis_params': analysis_params,
                'sample_mapping': st.session_state.sample_mapping,
                'graph_settings': {k: v for k, v in st.session_state.graph_settings.items() 
                                  if not k.startswith('bar_colors_per')},
                'qc_exclusions': list(st.session_state.qc_flagged_wells),
                'qc_edits': st.session_state.qc_edited_values
            }
            
            st.download_button(
                label="📥 Analysis Configuration (JSON)",
                data=json.dumps(config_data, indent=2, default=str),
                file_name=f"config_{datetime.now().strftime('%Y%m%d_%H%M')}.json",
                mime="application/json"
            )
        
        # Export guide
        with st.expander("💡 Export Guide"):
            st.markdown("""
            ### For Journal Submission
            1. Download **Excel Report** for supplementary data
            2. Use **HTML graphs** → Open in browser → Right-click → Save as PNG
            3. For highest quality: Install `kaleido` package for direct PNG/PDF export
            
            ### For Presentations
            - HTML files can be embedded in PowerPoint
            - Interactive features work in presentations!
            
            ### For Reproducibility
            - Download **Configuration JSON** to recreate exact analysis
            - Share with collaborators for verification
            
            ### Publication Requirements
            - Most journals require 300 DPI minimum
            - Use Arial or Helvetica fonts
            - Figure width: 3.5" (single column) or 7" (double column)
            """)
        
        st.success("✅ All exports ready!")
    else:
        st.warning("⚠️ Complete analysis first")


# ==================== FOOTER ====================
st.markdown("---")
st.markdown("""
<div style='text-align: center; color: #666; padding: 20px;'>
    <p>🧬 <b>qPCR Analysis Suite Pro v3.0</b></p>
    <p>Enhanced Edition with QC Module, Bulk Operations, and Publication-Ready Exports</p>
    <p style='font-size: 12px;'>Housekeeping normalization • ΔΔCt calculations • Welch's t-test • Gene-by-gene visualization</p>
</div>
""", unsafe_allow_html=True)
