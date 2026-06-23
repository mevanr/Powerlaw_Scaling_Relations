  def create_visualizations(self, plot_data, performance, ranking):
    """Create benchmark plots matching R style"""
    
    import matplotlib.pyplot as plt
    import seaborn as sns
    import numpy as np
    import os
    
    # Set style to match R's default
    plt.style.use('default')
    sns.set_style("whitegrid")
    
    # Create output directory
    os.makedirs(self.config.out_dir, exist_ok=True)
    
    # ============================================================================
    # FIGURE 1: Boxplots (matching Figure1_boxplots_with_keras.png)
    # ============================================================================
    
    # Filter methods to match R figure
    r_methods = ['ols', 'theilsen', 'siegel', 'nlme', 'majoraxis', 'rma', 
                 'deming_std', 'deming_wtd', 'ensemble', 'lmer', 'lmer_slopes']
    
    fig, axes = plt.subplots(1, 3, figsize=(15, 6))
    
    colors = ['#E41A1C', '#377EB8', '#4DAF4A']  # R-style colors
    
    for idx, beta in enumerate(self.config.beta_values):
        data = plot_data[
            (plot_data['beta_true'] == beta) & 
            (plot_data['method'].isin(r_methods))
        ].copy()
        
        # Reorder methods to match R figure
        method_order = ['ols', 'theilsen', 'siegel', 'nlme', 'majoraxis', 
                       'rma', 'deming_std', 'deming_wtd', 'ensemble', 
                       'lmer', 'lmer_slopes']
        data['method'] = pd.Categorical(data['method'], categories=method_order, ordered=True)
        data = data.sort_values('method')
        
        # Create boxplot
        bp = axes[idx].boxplot(
            [data[data['method'] == m]['estimate'].values for m in method_order if m in data['method'].unique()],
            patch_artist=True,
            boxprops=dict(facecolor=colors[idx], color='black', alpha=0.7),
            medianprops=dict(color='black', linewidth=2),
            whiskerprops=dict(color='black'),
            capprops=dict(color='black'),
            flierprops=dict(marker='o', markerfacecolor='gray', markersize=4, alpha=0.5)
        )
        
        # Add true value line
        axes[idx].axhline(y=beta, color='red', linestyle='--', linewidth=2, label=f'True β={beta}')
        
        # Customize
        axes[idx].set_title(f'True β = {beta}', fontsize=12, fontweight='bold')
        axes[idx].set_xlabel('Method', fontsize=10)
        axes[idx].set_ylabel('Estimated Exponent (β̂)', fontsize=10)
        axes[idx].set_xticklabels(method_order, rotation=45, ha='right')
        axes[idx].grid(True, alpha=0.3)
        axes[idx].legend(loc='upper right')
        
        # Set y-axis limits to match R figure
        if beta == 0.6:
            axes[idx].set_ylim(0.35, 0.85)
        elif beta == 1.0:
            axes[idx].set_ylim(0.75, 1.25)
        else:  # beta == 1.4
            axes[idx].set_ylim(1.15, 1.65)
    
    plt.suptitle('Figure 1: Comparison of Power-Law Exponent Estimation Methods\nRandom slopes (SD=0.40), Habitat-specific X means, σx=0.20, σy=0.40', 
                 fontsize=14, fontweight='bold')
    plt.tight_layout()
    plt.savefig(os.path.join(self.config.out_dir, 'Figure1_boxplots.png'), 
                dpi=300, bbox_inches='tight')
    plt.show()
    
    # ============================================================================
    # FIGURE 2: Bias Bar Plot (matching Figure2_bias_with_keras.png)
    # ============================================================================
    
    fig, ax = plt.subplots(figsize=(14, 6))
    
    # Filter and prepare bias data
    bias_methods = ['deming_std', 'deming_wtd', 'ensemble', 'lmer', 'lmer_slopes', 
                    'majoraxis', 'nlme', 'ols', 'rma', 'siegel', 'theilsen']
    
    bias_data = performance[performance['method'].isin(bias_methods)].copy()
    
    # Pivot for plotting
    pivot_bias = bias_data.pivot(index='method', columns='beta_true', values='bias')
    pivot_bias = pivot_bias.reindex(bias_methods)  # Reorder
    
    # Create grouped bar plot
    x = np.arange(len(pivot_bias.index))
    width = 0.25
    
    colors = ['#E41A1C', '#377EB8', '#4DAF4A']
    
    for i, beta in enumerate(self.config.beta_values):
        bars = ax.bar(x + i*width, pivot_bias[beta], width, 
                      label=f'True β={beta}', color=colors[i], alpha=0.8)
        
        # Add error bars (standard error)
        se_data = bias_data[bias_data['beta_true'] == beta].set_index('method')['std_est']
        ax.errorbar(x + i*width, pivot_bias[beta], yerr=se_data[pivot_bias.index]/np.sqrt(100),
                   fmt='none', color='black', capsize=3)
    
    ax.axhline(y=0, color='black', linestyle='--', linewidth=1)
    ax.set_xlabel('Method', fontsize=12)
    ax.set_ylabel('Mean Bias (β̂ - β)', fontsize=12)
    ax.set_title('Figure 2: Bias by Method (Random Slopes + Habitat-Specific X)', 
                 fontsize=14, fontweight='bold')
    ax.set_xticks(x + width)
    ax.set_xticklabels(pivot_bias.index, rotation=45, ha='right')
    ax.legend(loc='upper right')
    ax.grid(True, alpha=0.3, axis='y')
    
    plt.tight_layout()
    plt.savefig(os.path.join(self.config.out_dir, 'Figure2_bias.png'), 
                dpi=300, bbox_inches='tight')
    plt.show()
    
    # ============================================================================
    # FIGURE 3: RMSE Heatmap (matching Figure3_rmse_heatmap_with_keras.png)
    # ============================================================================
    
    fig, ax = plt.subplots(figsize=(12, 4))
    
    # Filter methods
    rmse_methods = ['deming_std', 'deming_wtd', 'ensemble', 'lmer', 'lmer_slopes', 
                    'majoraxis', 'nlme', 'ols', 'rma', 'theilsen', 'siegel']
    
    rmse_data = performance[performance['method'].isin(rmse_methods)].copy()
    
    # Pivot for heatmap
    pivot_rmse = rmse_data.pivot(index='method', columns='beta_true', values='rmse')
    pivot_rmse = pivot_rmse.reindex(rmse_methods)
    
    # Create heatmap
    im = ax.imshow(pivot_rmse.values, cmap='YlOrRd', aspect='auto', vmin=0.2, vmax=1.0)
    
    # Add text annotations
    for i in range(len(pivot_rmse.index)):
        for j in range(len(pivot_rmse.columns)):
            text = ax.text(j, i, f'{pivot_rmse.values[i, j]:.3f}',
                          ha='center', va='center', color='black' if pivot_rmse.values[i, j] < 0.6 else 'white')
    
    # Customize
    ax.set_xticks(range(len(pivot_rmse.columns)))
    ax.set_xticklabels([f'True β={b}' for b in pivot_rmse.columns])
    ax.set_yticks(range(len(pivot_rmse.index)))
    ax.set_yticklabels(pivot_rmse.index)
    ax.set_title('Figure 3: RMSE Heatmap (Random Slopes + Habitat-Specific X)', 
                 fontsize=14, fontweight='bold')
    
    plt.colorbar(im, ax=ax, label='RMSE')
    plt.tight_layout()
    plt.savefig(os.path.join(self.config.out_dir, 'Figure3_rmse_heatmap.png'), 
                dpi=300, bbox_inches='tight')
    plt.show()
    
    # ============================================================================
    # FIGURE 4: Method Ranking (matching Figure4_ranking_with_keras.png)
    # ============================================================================
    
    fig, ax = plt.subplots(figsize=(10, 6))
    
    # Use ranking data
    ranking_df = ranking.reset_index()
    ranking_df.columns = ['method', 'avg_rmse']
    
    # Filter to match R figure
    rank_methods = ['deming_std', 'deming_wtd', 'lmer', 'lmer_slopes', 'nlme', 
                    'majoraxis', 'theilsen', 'siegel', 'ols', 'ensemble', 'rma']
    
    ranking_df = ranking_df[ranking_df['method'].isin(rank_methods)]
    ranking_df = ranking_df.sort_values('avg_rmse', ascending=True)
    
    # Create horizontal bar plot
    colors = plt.cm.viridis(np.linspace(0, 1, len(ranking_df)))
    
    bars = ax.barh(range(len(ranking_df)), ranking_df['avg_rmse'], color=colors)
    
    # Add rank labels
    for i, (bar, (_, row)) in enumerate(zip(bars, ranking_df.iterrows())):
        ax.text(bar.get_width() + 0.01, bar.get_y() + bar.get_height()/2, 
                f'Rank {i+1}\n{row["avg_rmse"]:.3f}', 
                va='center', fontsize=9)
    
    ax.set_yticks(range(len(ranking_df)))
    ax.set_yticklabels(ranking_df['method'])
    ax.set_xlabel('Average RMSE', fontsize=12)
    ax.set_title('Figure 4: Method Ranking by Average RMSE\nWith random slopes, habitat-specific X means', 
                 fontsize=14, fontweight='bold')
    ax.invert_yaxis()  # Best at top
    ax.grid(True, alpha=0.3, axis='x')
    
    plt.tight_layout()
    plt.savefig(os.path.join(self.config.out_dir, 'Figure4_ranking.png'), 
                dpi=300, bbox_inches='tight')
    plt.show()
    
    return fig