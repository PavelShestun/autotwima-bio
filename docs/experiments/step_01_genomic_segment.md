# Эксперимент 1.1: Модель геномного сегмента в дискретно-бесконечной Phi-группе

## 1. Теоретическая основа
- **Модель статьи:** Раздел 5.2, $M_\infty = \mathbb{Z} \times \text{Bool}$.
- **Биоинформатическая редукция:** 
  - Геномная координата/длина $x \in \mathbb{Z}_{\ge 0}$.
  - Ориентация цепи $u \in \{F, T\}$, где Forward = $F$ ($\text{sgn} = +1$), Reverse = $T$ ($\text{sgn} = -1$).

## 2. Что реализовано
- Датакласс `GenomicSegment` с валидацией неотрицательной длины.
- Метод `to_phi()` детерминированного отображения в `PhiInfElement`.

## 3. Результаты тестов
- `test_genomic_segment_forward_strand_mapping`: passed
- `test_genomic_segment_reverse_strand_mapping`: passed
- `test_genomic_segment_rejects_negative_length`: passed
- **Quality Gates:** Ruff 0 ошибок, ty 100% type check passed, pytest 11/11 passed.
