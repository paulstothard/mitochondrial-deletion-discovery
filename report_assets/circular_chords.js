function supportSliderValueForThreshold(target, supportMin, supportMax, sliderMax) {
  const maxValue = Number(sliderMax);
  if (!Number.isFinite(target) || !Number.isFinite(supportMin) || !Number.isFinite(supportMax)
    || !Number.isFinite(maxValue) || maxValue <= 0 || supportMin <= 0) return 0;
  if (supportMax <= supportMin || target <= supportMin) return 0;
  if (target >= supportMax) return maxValue;
  const fraction = Math.log(target / supportMin) / Math.log(supportMax / supportMin);
  return Math.max(0, Math.min(maxValue, Math.floor(fraction * maxValue)));
}

function syncSupportSliderToObservations(entries, minimumObservations, minimumSize, slider,
  observationField, supportMin, supportMax) {
  const candidates = entries
    .filter((entry) => Number(entry.dataset[observationField]) >= minimumObservations
      && Number(entry.dataset.deletedSize) >= minimumSize)
    .map((entry) => Number(entry.dataset.support))
    .filter(Number.isFinite);
  const target = candidates.length ? Math.min(...candidates) : supportMax;
  slider.value = String(candidates.length
    ? supportSliderValueForThreshold(target, supportMin, supportMax, slider.max)
    : Number(slider.max));
}

document.querySelectorAll('[data-rainfall-controls]').forEach((controls) => {
  const target = document.getElementById(controls.dataset.target);
  const points = Array.from(target.querySelectorAll('.rainfall-point'));
  const slider = controls.querySelector('[data-rainfall-support-slider]');
  const supportOutput = controls.querySelector('[data-rainfall-support-output]');
  const observationFilter = controls.querySelector('[data-observation-filter]');
  const linkedOption = controls.querySelector('[data-linked-option]');
  const sizeFilter = controls.querySelector('[data-size-filter]');
  const reset = controls.querySelector('[data-reset-rainfall-controls]');
  const status = controls.querySelector('[data-rainfall-filter-status]');
  const supports = points.map((point) => Number(point.dataset.support)).filter(Number.isFinite);
  const supportMin = Math.min(...supports);
  const supportMax = Math.max(...supports);

  function filteredObservationMinimum(pointsToCheck) {
    const observations = pointsToCheck
      .map((point) => Number(point.dataset.supportingReads))
      .filter(Number.isFinite);
    return observations.length ? Math.min(...observations) : 0;
  }

  function formatRainfallSupport(value) {
    if (value >= 100) return value.toFixed(0);
    if (value >= 10) return value.toFixed(1).replace(/\.0$/, '');
    if (value >= 1) return value.toFixed(2).replace(/0+$/, '').replace(/\.$/, '');
    return value.toPrecision(3).replace(/0+$/, '').replace(/\.$/, '');
  }

  function supportThreshold() {
    const fraction = Number(slider.value) / Number(slider.max);
    if (!Number.isFinite(supportMin) || !Number.isFinite(supportMax) || supportMin <= 0) return 0;
    if (fraction <= 0) return 0;
    if (supportMax <= supportMin) return supportMin;
    return supportMin * Math.pow(supportMax / supportMin, fraction);
  }

  function renderRainfall() {
    const threshold = supportThreshold();
    const minimumSize = sizeFilter ? Number(sizeFilter.value) : 0;
    const supportAndSizeEligible = points.filter((point) => Number(point.dataset.support) >= threshold
      && Number(point.dataset.deletedSize) >= minimumSize);
    const linkedMinimum = filteredObservationMinimum(supportAndSizeEligible);
    const minimumObservations = observationFilter && observationFilter.value !== 'linked'
      ? Number(observationFilter.value)
      : linkedMinimum;
    let visible = 0;
    points.forEach((point) => {
      const keep = Number(point.dataset.support) >= threshold
        && Number(point.dataset.deletedSize) >= minimumSize
        && Number(point.dataset.supportingReads) >= minimumObservations;
      point.style.display = keep ? '' : 'none';
      if (keep) visible += 1;
    });
    if (linkedOption) {
      linkedOption.textContent = linkedMinimum > 0
        ? `Auto (\u2265 ${linkedMinimum.toLocaleString()})`
        : 'Auto (no calls pass)';
    }
    supportOutput.textContent = Number(slider.value) === 0
      ? 'All loaded calls'
      : `>= ${formatRainfallSupport(threshold)}`;
    const sizeNote = minimumSize > 0 ? `, deleted size >= ${minimumSize.toLocaleString()} bp` : '';
    const observationNote = minimumObservations > 0 ? `, supporting observations >= ${minimumObservations.toLocaleString()}` : '';
    status.textContent = `Showing ${visible.toLocaleString()} of ${points.length.toLocaleString()} loaded exact deletions${sizeNote}${observationNote}`;
  }

  slider.addEventListener('input', () => {
    if (observationFilter) observationFilter.value = 'linked';
    renderRainfall();
  });
  function syncRainfallSupportToObservations() {
    if (!observationFilter || observationFilter.value === 'linked') return;
    const minimumSize = sizeFilter ? Number(sizeFilter.value) : 0;
    syncSupportSliderToObservations(
      points, Number(observationFilter.value), minimumSize, slider,
      'supportingReads', supportMin, supportMax,
    );
  }

  if (observationFilter) {
    observationFilter.addEventListener('change', () => {
      syncRainfallSupportToObservations();
      renderRainfall();
    });
  }
  if (sizeFilter) {
    sizeFilter.addEventListener('change', () => {
      syncRainfallSupportToObservations();
      renderRainfall();
    });
  }
  reset.addEventListener('click', () => {
    slider.value = '0';
    if (observationFilter) observationFilter.value = 'linked';
    if (sizeFilter) sizeFilter.value = '0';
    renderRainfall();
  });
  renderRainfall();
});


    document.querySelectorAll('[data-chord-controls]').forEach((controls) => {
      const target = document.getElementById(controls.dataset.target);
      const chords = Array.from(target.querySelectorAll('.deletion-chord'));
      const slider = controls.querySelector('[data-support-slider]');
      const supportOutput = controls.querySelector('[data-support-output]');
      const observationFilter = controls.querySelector('[data-observation-filter]');
      const sizeFilter = controls.querySelector('[data-size-filter]');
      const linkedOption = controls.querySelector('[data-linked-option]');
      const reset = controls.querySelector('[data-reset-controls]');
      const status = controls.querySelector('[data-filter-status]');
      const supports = chords.map((chord) => Number(chord.dataset.support)).filter(Number.isFinite);
      const supportMin = Math.min(...supports);
      const supportMax = Math.max(...supports);
      let baselineMode = false;

      function formatSupport(value) {
        if (value >= 100) return value.toFixed(0);
        if (value >= 10) return value.toFixed(1).replace(/\.0$/, '');
        if (value >= 1) return value.toFixed(2).replace(/0+$/, '').replace(/\.$/, '');
        return value.toPrecision(3).replace(/0+$/, '').replace(/\.$/, '');
      }

      function supportThreshold() {
        const fraction = Number(slider.value) / Number(slider.max);
        if (!Number.isFinite(supportMin) || !Number.isFinite(supportMax) || supportMin <= 0) return 0;
        if (fraction <= 0 || supportMax <= supportMin) return supportMin;
        return supportMin * Math.pow(supportMax / supportMin, fraction);
      }

      function syncChordSupportToObservations() {
        if (observationFilter.value === 'linked') return;
        const minimumSize = sizeFilter ? Number(sizeFilter.value) : 0;
        syncSupportSliderToObservations(
          chords, Number(observationFilter.value), minimumSize, slider,
          'observations', supportMin, supportMax,
        );
      }

      function render() {
        const threshold = supportThreshold();
        const minimumSize = sizeFilter ? Number(sizeFilter.value) : 0;
        const linked = observationFilter.value === 'linked';
        const supportEligible = chords.filter((chord) => Number(chord.dataset.support) >= threshold);
        const sizeEligible = supportEligible.filter((chord) => Number(chord.dataset.deletedSize) >= minimumSize);
        const candidateChords = baselineMode
          ? chords.filter((chord) => chord.dataset.baseline === '1')
          : sizeEligible;
        const candidateObservations = candidateChords.map((chord) => Number(chord.dataset.observations));
        const linkedMinimum = candidateObservations.length ? Math.min(...candidateObservations) : 0;
        const minimumObservations = linked ? linkedMinimum : Number(observationFilter.value);
        let visible = 0;
        chords.forEach((chord) => {
          const passesPrimary = baselineMode
            ? chord.dataset.baseline === '1'
            : Number(chord.dataset.support) >= threshold
              && Number(chord.dataset.deletedSize) >= minimumSize;
          const keep = passesPrimary && Number(chord.dataset.observations) >= minimumObservations;
          chord.style.display = keep ? '' : 'none';
          if (keep) visible += 1;
        });
        linkedOption.textContent = linkedMinimum > 0
          ? `Auto (\u2265 ${linkedMinimum.toLocaleString()})`
          : 'Auto (no calls pass)';
        supportOutput.textContent = baselineMode
          ? 'PDF baseline'
          : (Number(slider.value) === 0 ? 'All loaded calls' : `>= ${formatSupport(threshold)}`);
        const baselineNote = baselineMode ? ' (PDF baseline)' : '';
        const sizeNote = minimumSize > 0 ? `, deleted size >= ${minimumSize.toLocaleString()} bp` : '';
        status.textContent = `Showing ${visible.toLocaleString()} of ${chords.length.toLocaleString()} loaded exact deletions${sizeNote}${baselineNote}`;
      }

      slider.addEventListener('input', () => {
        baselineMode = false;
        observationFilter.value = 'linked';
        render();
      });
      observationFilter.addEventListener('change', () => {
        baselineMode = false;
        syncChordSupportToObservations();
        render();
      });
      if (sizeFilter) {
        sizeFilter.addEventListener('change', () => {
          baselineMode = false;
          syncChordSupportToObservations();
          render();
        });
      }
      reset.addEventListener('click', () => {
        baselineMode = true;
        slider.value = '0';
        observationFilter.value = 'linked';
        if (sizeFilter) sizeFilter.value = '0';
        render();
      });
      render();
    });

    document.querySelectorAll('[data-breakpoint-pair-controls]').forEach((controls) => {
      const target = document.getElementById(controls.dataset.target);
      const points = Array.from(target.querySelectorAll('.breakpoint-pair-point'));
      const rankLabels = Array.from(target.querySelectorAll('[id^="breakpoint-pair-rank-"]'));
      const slider = controls.querySelector('[data-support-slider]');
      const supportOutput = controls.querySelector('[data-support-output]');
      const observationFilter = controls.querySelector('[data-observation-filter]');
      const linkedOption = controls.querySelector('[data-linked-option]');
      const sizeFilter = controls.querySelector('[data-size-filter]');
      const reset = controls.querySelector('[data-reset-breakpoint-pair-controls]');
      const status = controls.querySelector('[data-filter-status]');
      const supports = points.map((point) => Number(point.dataset.support)).filter(Number.isFinite);
      const supportMin = Math.min(...supports);
      const supportMax = Math.max(...supports);

      function formatSupport(value) {
        if (value >= 100) return value.toFixed(0);
        if (value >= 10) return value.toFixed(1).replace(/\.0$/, '');
        if (value >= 1) return value.toFixed(2).replace(/0+$/, '').replace(/\.$/, '');
        return value.toPrecision(3).replace(/0+$/, '').replace(/\.$/, '');
      }

      function supportThreshold() {
        const fraction = Number(slider.value) / Number(slider.max);
        if (!Number.isFinite(supportMin) || !Number.isFinite(supportMax) || supportMin <= 0) return 0;
        if (fraction <= 0 || supportMax <= supportMin) return supportMin;
        return supportMin * Math.pow(supportMax / supportMin, fraction);
      }

      function syncBreakpointSupportToObservations() {
        if (observationFilter.value === 'linked') return;
        const minimumSize = sizeFilter ? Number(sizeFilter.value) : 0;
        syncSupportSliderToObservations(
          points, Number(observationFilter.value), minimumSize, slider,
          'supportingObservations', supportMin, supportMax,
        );
      }

      function render() {
        const threshold = supportThreshold();
        const minimumSize = sizeFilter ? Number(sizeFilter.value) : 0;
        const supportAndSizeEligible = points.filter((point) => Number(point.dataset.support) >= threshold
          && Number(point.dataset.deletedSize) >= minimumSize);
        const eligibleObservations = supportAndSizeEligible
          .map((point) => Number(point.dataset.supportingObservations))
          .filter(Number.isFinite);
        const linkedMinimum = eligibleObservations.length ? Math.min(...eligibleObservations) : 0;
        const minimumObservations = observationFilter.value === 'linked'
          ? linkedMinimum
          : Number(observationFilter.value);
        const visibleRanks = new Set();
        let visible = 0;
        points.forEach((point) => {
          const keep = Number(point.dataset.support) >= threshold
            && Number(point.dataset.deletedSize) >= minimumSize
            && Number(point.dataset.supportingObservations) >= minimumObservations;
          point.style.display = keep ? '' : 'none';
          if (keep) {
            visible += 1;
            if (point.dataset.rank) visibleRanks.add(point.dataset.rank);
          }
        });
        rankLabels.forEach((label) => {
          const rank = label.id.replace('breakpoint-pair-rank-', '');
          label.style.display = visibleRanks.has(rank) ? '' : 'none';
        });
        linkedOption.textContent = linkedMinimum > 0
          ? `Auto (\u2265 ${linkedMinimum.toLocaleString()})`
          : 'Auto (no calls pass)';
        supportOutput.textContent = Number(slider.value) === 0
          ? 'All loaded calls'
          : `>= ${formatSupport(threshold)}`;
        const sizeNote = minimumSize > 0 ? `, deleted size >= ${minimumSize.toLocaleString()} bp` : '';
        const observationNote = minimumObservations > 0 ? `, supporting observations >= ${minimumObservations.toLocaleString()}` : '';
        status.textContent = `Showing ${visible.toLocaleString()} of ${points.length.toLocaleString()} loaded breakpoint pairs${sizeNote}${observationNote}`;
      }

      slider.addEventListener('input', () => {
        observationFilter.value = 'linked';
        render();
      });
      observationFilter.addEventListener('change', () => {
        syncBreakpointSupportToObservations();
        render();
      });
      if (sizeFilter) {
        sizeFilter.addEventListener('change', () => {
          syncBreakpointSupportToObservations();
          render();
        });
      }
      reset.addEventListener('click', () => {
        slider.value = '0';
        observationFilter.value = 'linked';
        if (sizeFilter) sizeFilter.value = '0';
        render();
      });
      render();
    });

    document.querySelectorAll('[data-comparison-controls]').forEach((controls) => {
      const target = document.getElementById(controls.dataset.target);
      const chords = Array.from(target.querySelectorAll('.comparison-chord'));
      const preset = controls.querySelector('[data-comparison-preset]');
      const presetGuidance = controls.querySelector('[data-comparison-preset-guidance]');
      const observationSlider = controls.querySelector('[data-comparison-observation-slider]');
      const observationOutput = controls.querySelector('[data-comparison-observation-output]');
      const differenceSlider = controls.querySelector('[data-comparison-difference-slider]');
      const differenceOutput = controls.querySelector('[data-comparison-difference-output]');
      const direction = controls.querySelector('[data-comparison-direction]');
      const resetRefinements = controls.querySelector('[data-reset-comparison-refinements]');
      const status = controls.querySelector('[data-comparison-status]');
      const observationValues = chords.map((chord) => Number(chord.dataset.totalObservations)).filter(Number.isFinite);
      const differenceValues = chords.map((chord) => Number(chord.dataset.absoluteDifference)).filter((value) => Number.isFinite(value) && value > 0);
      const presetRules = {
        'all': {field: null, threshold: 1},
        'replicate-significant': {field: 'replicateQ', threshold: 0.05},
        'replicate-suggestive': {field: 'replicateP', threshold: 0.05},
        'depth-significant': {field: 'depthQ', threshold: 0.05},
      };
      const observationMin = Math.max(1, Math.min(...observationValues));
      const observationMax = Math.max(observationMin, Math.max(...observationValues));
      const differenceMin = differenceValues.length ? Math.min(...differenceValues) : 0;
      const differenceMax = differenceValues.length ? Math.max(...differenceValues) : 0;

      function formatComparisonNumber(value) {
        if (!Number.isFinite(value)) return 'NA';
        if (value === 0) return '0';
        if (Math.abs(value) < 0.001) return value.toExponential(2);
        return value.toLocaleString(undefined, {maximumSignificantDigits: 3});
      }

      function optionalNumber(value) {
        if (value === undefined || value === null || String(value).trim() === '') return NaN;
        return Number(value);
      }

      function logThreshold(slider, minimum, maximum, zeroMeansAll = false) {
        const fraction = Number(slider.value) / Number(slider.max);
        if (zeroMeansAll && fraction <= 0) return 0;
        if (fraction >= 1) return maximum;
        if (!Number.isFinite(minimum) || !Number.isFinite(maximum) || minimum <= 0 || maximum <= minimum) return minimum;
        return minimum * Math.pow(maximum / minimum, fraction);
      }

      function renderComparison() {
        const rule = presetRules[preset.value] || presetRules.all;
        const observationThreshold = Math.ceil(logThreshold(observationSlider, observationMin, observationMax));
        const differenceThreshold = logThreshold(differenceSlider, differenceMin, differenceMax, true);
        let visible = 0;
        chords.forEach((chord) => {
          const statValue = rule.field === null ? 0 : optionalNumber(chord.dataset[rule.field]);
          const passesStatistic = rule.field === null || (Number.isFinite(statValue) && statValue <= rule.threshold);
          const change = Number(chord.dataset.difference);
          const passesDirection = direction.value === 'both'
            || (direction.value === 'right' && change > 0)
            || (direction.value === 'left' && change < 0);
          const keep = passesStatistic
            && Number(chord.dataset.totalObservations) >= observationThreshold
            && Number(chord.dataset.absoluteDifference) >= differenceThreshold
            && passesDirection;
          chord.style.display = keep ? '' : 'none';
          if (keep) visible += 1;
        });
        observationOutput.textContent = `\u2265 ${observationThreshold.toLocaleString()}`;
        differenceOutput.textContent = Number(differenceSlider.value) === 0
          ? 'All differences'
          : `\u2265 ${formatComparisonNumber(differenceThreshold)}`;
        status.textContent = `Showing ${visible.toLocaleString()} of ${chords.length.toLocaleString()} exact deletion comparisons`;
      }

      function applyPreset() {
        if (preset.value === 'replicate-significant') {
          presetGuidance.textContent = 'Use this view for biological group conclusions. BH q \u2264 0.05 accounts for testing many exact deletions; an empty plot means none pass this threshold.';
        } else if (preset.value === 'replicate-suggestive') {
          presetGuidance.textContent = 'Shows unadjusted replicate-level p \u2264 0.05 results for exploration. These are not significant after correction unless they also pass the BH-q view.';
        } else if (preset.value === 'depth-significant') {
          presetGuidance.textContent = 'Shows BH-adjusted read-depth enrichment. This is technical read-count evidence and must not be described as biological-replicate significance.';
        } else {
          presetGuidance.textContent = 'Shows every delivered exact-deletion comparison. Use Replicate-significant for biological group conclusions.';
        }
        renderComparison();
      }

      [observationSlider, differenceSlider].forEach((slider) => slider.addEventListener('input', renderComparison));
      direction.addEventListener('change', renderComparison);
      preset.addEventListener('change', applyPreset);
      resetRefinements.addEventListener('click', () => {
        observationSlider.value = '0';
        differenceSlider.value = '0';
        direction.value = 'both';
        renderComparison();
      });
      applyPreset();
    });

    const hoverTooltip = document.createElement('div');
    hoverTooltip.className = 'hover-tooltip';
    hoverTooltip.setAttribute('role', 'tooltip');
    document.body.appendChild(hoverTooltip);

    function addTooltipRow(label, value) {
      const row = document.createElement('div');
      row.className = 'hover-tooltip-row';
      const prefix = document.createElement('span');
      prefix.textContent = `${label}: `;
      row.append(prefix, document.createTextNode(value));
      hoverTooltip.appendChild(row);
    }

    function formatTooltipNumber(value) {
      if (value === undefined || value === null || String(value).trim() === '') return 'NA';
      const number = Number(value);
      if (!Number.isFinite(number)) return 'NA';
      if (number !== 0 && Math.abs(number) < 0.001) return number.toExponential(3);
      return number.toLocaleString(undefined, {maximumSignificantDigits: 4});
    }

    function formatTooltipValue(value) {
      if (value === undefined || value === null || String(value).trim() === '') return 'NA';
      const number = Number(value);
      return Number.isFinite(number) ? formatTooltipNumber(number) : String(value);
    }

    function populateTooltip(target) {
      hoverTooltip.replaceChildren();
      const heading = document.createElement('strong');
      if (target.classList.contains('rainfall-point')) {
        heading.textContent = target.dataset.exactDeletionId || 'Exact deletion';
        hoverTooltip.appendChild(heading);
        addTooltipRow('Group', target.dataset.group || 'NA');
        addTooltipRow('Directed breakpoints', `${formatTooltipNumber(target.dataset.leftBreakpoint)} to ${formatTooltipNumber(target.dataset.rightBreakpoint)}`);
        addTooltipRow('Deleted size', `${formatTooltipNumber(target.dataset.deletedSize)} bp`);
        addTooltipRow(target.dataset.supportLabel || 'Plotted support', formatTooltipNumber(target.dataset.support));
        addTooltipRow('Supporting observations', formatTooltipNumber(target.dataset.supportingReads));
        addTooltipRow('Affected features', (target.dataset.affectedFeatures || 'NA').replaceAll('_', ' '));
        addTooltipRow('Arc annotation', (target.dataset.arcContext || 'NA').replaceAll('_', ' '));
        addTooltipRow('Major/minor arc bp', `${formatTooltipNumber(target.dataset.majorArcBp)} / ${formatTooltipNumber(target.dataset.minorArcBp)}`);
        addTooltipRow('Origin-spanning', target.dataset.crossesOrigin || 'NA');
        if (target.dataset.knownDeletion) addTooltipRow('Configured match', target.dataset.knownDeletion.replaceAll('_', ' '));
      } else if (target.classList.contains('endpoint-density-bin')) {
        heading.textContent = `${formatTooltipNumber(target.dataset.binStart)}-${formatTooltipNumber(target.dataset.binEnd)} bp`;
        hoverTooltip.appendChild(heading);
        addTooltipRow('Group', target.dataset.group || 'NA');
        addTooltipRow('Bin midpoint', `${formatTooltipNumber(target.dataset.binMidpoint)} bp`);
        const supportLabel = target.dataset.supportLabel || 'Plotted support';
        addTooltipRow('Left exact deletion calls', formatTooltipNumber(target.dataset.leftEndpointCount));
        addTooltipRow(`Left ${supportLabel}`, formatTooltipNumber(target.dataset.leftSupport));
        addTooltipRow('Left raw supporting observations', formatTooltipNumber(target.dataset.leftRawSupportingReads));
        addTooltipRow('Right exact deletion calls', formatTooltipNumber(target.dataset.rightEndpointCount));
        addTooltipRow(`Right ${supportLabel}`, formatTooltipNumber(target.dataset.rightSupport));
        addTooltipRow('Right raw supporting observations', formatTooltipNumber(target.dataset.rightRawSupportingReads));
        addTooltipRow('Total distinct exact deletion calls', formatTooltipNumber(target.dataset.endpointCount));
        addTooltipRow('Total raw supporting observations', formatTooltipNumber(target.dataset.rawSupportingReads));
        addTooltipRow(`Total ${supportLabel}`, formatTooltipNumber(target.dataset.summedSupport));
        addTooltipRow(`Smoothed ${supportLabel}`, formatTooltipNumber(target.dataset.smoothedSupport));
        addTooltipRow('Smoothed exact deletion call count', formatTooltipNumber(target.dataset.smoothedEndpointCount));
      } else if (target.classList.contains('breakpoint-pair-point')) {
        heading.textContent = target.dataset.exactDeletionId || 'Breakpoint pair';
        hoverTooltip.appendChild(heading);
        addTooltipRow('Group', target.dataset.group || 'NA');
        addTooltipRow('Directed breakpoints', `${formatTooltipNumber(target.dataset.leftBreakpoint)} to ${formatTooltipNumber(target.dataset.rightBreakpoint)}`);
        addTooltipRow('Deleted size', `${formatTooltipNumber(target.dataset.deletedSize)} bp`);
        addTooltipRow(target.dataset.supportLabel || 'Plotted support', formatTooltipNumber(target.dataset.support));
        addTooltipRow('Supporting observations', formatTooltipNumber(target.dataset.supportingObservations));
        addTooltipRow('Exact deletions in pair', formatTooltipNumber(target.dataset.pairCount));
        addTooltipRow('Affected features', (target.dataset.affectedFeatures || 'NA').replaceAll('_', ' '));
        addTooltipRow('Origin-spanning', target.dataset.crossesOrigin || 'NA');
        addTooltipRow('Arc annotation', (target.dataset.arcContext || 'NA').replaceAll('_', ' '));
      } else if (target.classList.contains('group-mean-point')) {
        heading.textContent = `Group mean: ${target.dataset.group || 'NA'}`;
        hoverTooltip.appendChild(heading);
        addTooltipRow('Group', target.dataset.group || 'NA');
        addTooltipRow(target.dataset.yLabel || 'Mean value', formatTooltipValue(target.dataset.yValue));
        addTooltipRow('Samples contributing', formatTooltipNumber(target.dataset.sampleCount));
        addTooltipRow('Sample IDs', target.dataset.samples || 'NA');
        if (target.dataset.ciLow && target.dataset.ciHigh) {
          addTooltipRow('95% CI', `${formatTooltipValue(target.dataset.ciLow)} to ${formatTooltipValue(target.dataset.ciHigh)}`);
        }
        if (target.dataset.age) addTooltipRow('Age', target.dataset.age);
        if (target.dataset.treatment) addTooltipRow('Treatment', target.dataset.treatment);
      } else if (target.classList.contains('bar-plot-bar')) {
        const binStart = target.dataset.binStart;
        const binEnd = target.dataset.binEnd;
        heading.textContent = binStart && binEnd
          ? `${formatTooltipValue(binStart)}-${formatTooltipValue(binEnd)} bp`
          : (target.dataset.category || target.dataset.feature || target.dataset.label || 'Bar');
        hoverTooltip.appendChild(heading);
        if (target.dataset.group) addTooltipRow('Group', target.dataset.group);
        if (target.dataset.category) addTooltipRow('Category', target.dataset.category.replaceAll('_', ' '));
        if (target.dataset.feature) addTooltipRow('Feature', target.dataset.feature.replaceAll('_', ' '));
        if (target.dataset.label && target.dataset.label !== target.dataset.category) addTooltipRow('Displayed label', target.dataset.label);
        if (binStart && binEnd) addTooltipRow('Size interval', `${formatTooltipValue(binStart)} to ${formatTooltipValue(binEnd)} bp`);
        if (target.dataset.deletionId) addTooltipRow('Exact deletion', target.dataset.deletionId);
        if (target.dataset.groupValues) {
          try {
            const groupValues = JSON.parse(target.dataset.groupValues);
            groupValues.forEach((entry) => {
              const group = entry.group || 'Group';
              addTooltipRow(`${group}: ${entry.valueLabel || 'Plotted value'}`, formatTooltipValue(entry.value));
              const supportingReads = entry.supportingReads ?? entry.supporting_reads;
              addTooltipRow(`${group}: supporting reads`, formatTooltipNumber(supportingReads));
            });
          } catch (error) {
            addTooltipRow('Group values', target.dataset.groupValues);
          }
        } else {
          if (target.dataset.supportingReads) addTooltipRow('Supporting reads', formatTooltipNumber(target.dataset.supportingReads));
          if (target.dataset.value) addTooltipRow(target.dataset.valueLabel || 'Plotted value', formatTooltipValue(target.dataset.value));
        }
      } else if (target.classList.contains('ordination-point') || target.classList.contains('sample-point')) {
        heading.textContent = target.dataset.sample || 'Sample';
        hoverTooltip.appendChild(heading);
        addTooltipRow('Group', target.dataset.group || 'NA');
        addTooltipRow(target.dataset.xLabel || 'X', formatTooltipValue(target.dataset.xValue));
        addTooltipRow(target.dataset.yLabel || 'Y', formatTooltipValue(target.dataset.yValue));
        if (target.dataset.biologicalReplicate) addTooltipRow('Biological replicate', target.dataset.biologicalReplicate);
        if (target.dataset.layout) addTooltipRow('Read layout', target.dataset.layout);
        if (target.dataset.tissue) addTooltipRow('Tissue', target.dataset.tissue);
        if (target.dataset.age) addTooltipRow('Age', target.dataset.age);
        if (target.dataset.treatment) addTooltipRow('Treatment', target.dataset.treatment);
      } else if (target.classList.contains('comparison-chord')) {
        heading.textContent = `Comparison rank ${target.dataset.rank}: ${target.dataset.deletionId}`;
        hoverTooltip.appendChild(heading);
        addTooltipRow('Directed breakpoints', `${Number(target.dataset.leftBreakpoint).toLocaleString()} to ${Number(target.dataset.rightBreakpoint).toLocaleString()}`);
        addTooltipRow('Deleted interval', `${Number(target.dataset.deletedSize).toLocaleString()} bp`);
        addTooltipRow(`${target.dataset.leftGroup} mean`, formatTooltipNumber(target.dataset.leftMean));
        addTooltipRow(`${target.dataset.rightGroup} mean`, formatTooltipNumber(target.dataset.rightMean));
        addTooltipRow(`${target.dataset.rightGroup} minus ${target.dataset.leftGroup}`, formatTooltipNumber(target.dataset.difference));
        addTooltipRow('Supporting observations', `${Number(target.dataset.leftObservations).toLocaleString()} / ${Number(target.dataset.rightObservations).toLocaleString()}`);
        addTooltipRow('Samples with signal', formatTooltipNumber(target.dataset.samplesWithSignal));
        addTooltipRow('Replicate p / BH q', `${formatTooltipNumber(target.dataset.replicateP)} / ${formatTooltipNumber(target.dataset.replicateQ)}`);
        addTooltipRow('Read-depth Fisher p / BH q', `${formatTooltipNumber(target.dataset.depthP)} / ${formatTooltipNumber(target.dataset.depthQ)}`);
        addTooltipRow('Affected features', (target.dataset.affectedFeatures || 'NA').replaceAll('_', ' '));
        addTooltipRow('Arc annotation', target.dataset.arcContext.replaceAll('_', ' '));
        if (target.dataset.knownDeletion) addTooltipRow('Known-deletion match', target.dataset.knownDeletion.replaceAll('_', ' '));
      } else if (target.classList.contains('deletion-chord')) {
        heading.textContent = `Rank ${target.dataset.rank}: ${target.dataset.deletionId}`;
        hoverTooltip.appendChild(heading);
        addTooltipRow('Directed breakpoints', `${Number(target.dataset.leftBreakpoint).toLocaleString()} to ${Number(target.dataset.rightBreakpoint).toLocaleString()}`);
        addTooltipRow('Deleted interval', `${Number(target.dataset.deletedSize).toLocaleString()} bp`);
        addTooltipRow(target.dataset.supportLabel || 'Plotted support', Number(target.dataset.support).toLocaleString(undefined, {maximumSignificantDigits: 4}));
        addTooltipRow('Supporting observations', Number(target.dataset.observations).toLocaleString());
        addTooltipRow('Affected features', (target.dataset.affectedFeatures || 'NA').replaceAll('_', ' '));
        addTooltipRow('Arc annotation', target.dataset.arcContext.replaceAll('_', ' '));
        addTooltipRow('Major/minor arc bp', `${Number(target.dataset.majorArcBp).toLocaleString()} / ${Number(target.dataset.minorArcBp).toLocaleString()}`);
      } else {
        heading.textContent = target.dataset.featureName;
        hoverTooltip.appendChild(heading);
        addTooltipRow('Feature type', target.dataset.featureType);
        addTooltipRow('Coordinates', `${Number(target.dataset.featureStart).toLocaleString()} to ${Number(target.dataset.featureEnd).toLocaleString()}`);
      }
    }

    function positionTooltip(event) {
      const gap = 14;
      let left = event.clientX + gap;
      let top = event.clientY + gap;
      if (left + hoverTooltip.offsetWidth > window.innerWidth - 8) left = event.clientX - hoverTooltip.offsetWidth - gap;
      if (top + hoverTooltip.offsetHeight > window.innerHeight - 8) top = event.clientY - hoverTooltip.offsetHeight - gap;
      hoverTooltip.style.left = `${Math.max(8, left)}px`;
      hoverTooltip.style.top = `${Math.max(8, top)}px`;
    }

    document.addEventListener('pointerover', (event) => {
      const target = event.target instanceof Element
        ? event.target.closest('.rainfall-point, .endpoint-density-bin, .breakpoint-pair-point, .ordination-point, .sample-point, .group-mean-point, .bar-plot-bar, .deletion-chord, .comparison-chord, .mt-feature')
        : null;
      if (!target) return;
      populateTooltip(target);
      hoverTooltip.style.display = 'block';
      positionTooltip(event);
    });
    document.addEventListener('pointermove', (event) => {
      if (hoverTooltip.style.display === 'block') positionTooltip(event);
    });
    document.addEventListener('pointerout', (event) => {
      const target = event.target instanceof Element
        ? event.target.closest('.rainfall-point, .endpoint-density-bin, .breakpoint-pair-point, .ordination-point, .sample-point, .group-mean-point, .bar-plot-bar, .deletion-chord, .comparison-chord, .mt-feature')
        : null;
      const related = event.relatedTarget instanceof Element
        ? event.relatedTarget.closest('.rainfall-point, .endpoint-density-bin, .breakpoint-pair-point, .ordination-point, .sample-point, .group-mean-point, .bar-plot-bar, .deletion-chord, .comparison-chord, .mt-feature')
        : null;
      if (target && target !== related) hoverTooltip.style.display = 'none';
    });
