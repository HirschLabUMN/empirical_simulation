traits=c('EHT', 'PHT', 'Moisture', 'YLD')
n_markers=c(1, 3, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000)
effects=c(1, 0.1)
iters=c(1:3)
models=c('A', 'D')

complete_preds = tibble()

for(trait in traits){
  for(marker in n_markers){
    for(effect in effects){
      for(iter in iters){
        for(model in models){
          present = file.exists(paste0('~/Desktop/Grad_School/Research/empirical_sim/analysis/sim_traits/',
                                  trait,
                                  '/traits/avg_rank_all/n_markers_',
                                  marker,
                                  '/effects_',
                                  effect,
                                  '/prediction_results/iter',
                                  iter,
                                  '/model-',
                                  model,
                                  '/Gstr-diag_Rstr-us/prediction_accuracy.CV1.txt'))
          complete_preds = complete_preds %>%
            bind_rows(tibble(Trait = trait,
                             N_Marker = marker,
                             Effect = effect,
                             Iter = iter,
                             Model = model,
                             Present = present))
        }
      }
    }
  }
}

complete_preds %>%
  filter(Present == T)

complete_preds %>%
  group_by(Trait, N_Marker, Effect, Model) %>%
  summarise(Prop_Present = sum(Present) / n()) %>%
  ggplot(aes(x = factor(N_Marker), y = Trait, fill = Prop_Present))+
  geom_tile()+
  facet_grid(Effect ~ Model)+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

complete_preds %>%
  group_by(Trait, N_Marker, Effect, Model) %>%
  summarise(Prop_Present = sum(Present) / n()) %>%
  filter(Prop_Present == 0)


