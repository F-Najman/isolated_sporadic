-- Source columns for an audit on an LMFDB PostgreSQL mirror.
-- Bind the levels and orbit labels from the shipped data extracts.
SELECT label,level,dim,char_order,analytic_rank,analytic_rank_proved,
       CASE WHEN analytic_rank=0 AND analytic_rank_proved THEN NULL ELSE traces[1:20] END
  FROM mf_newforms
 WHERE weight=2 AND level=ANY(:levels)
 ORDER BY level,label;

SELECT label,level,hecke_ring_rank,hecke_ring_power_basis,
       hecke_ring_cyclotomic_generator,field_poly,hecke_ring_numerators,
       hecke_ring_denominators,hecke_ring_character_values,
       ap->0,ap->1,ap->2,ap->3,ap->4,ap->5,ap->6,ap->7
  FROM mf_hecke_nf WHERE label=ANY(:labels) ORDER BY label;
