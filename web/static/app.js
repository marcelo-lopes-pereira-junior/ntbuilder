/* ============================================================================
   NTBuilder — web frontend
   ----------------------------------------------------------------------------
   Vanilla JS state machine that drives the FastAPI backend.  The desktop
   feature set is mirrored here, including:
     - structure load (upload / example) and primitive-cell reduction
     - polar Hamada map (Plotly) with strain filter and atoms/strain coloring
     - single-tube build with 9 output formats
     - multi-wall (MWNT) with realised-wall-plan rendering
     - bundle and deform/torsion (chained via from_job_id)
     - structure analysis (bond histogram + electronic + symmetry)
     - ready-to-paste Methods text and DFT input generation
     - undo / redo over the job history (Ctrl+Z / Ctrl+Y)
     - chaining-warning prompts mirroring QMessageBox.question on desktop
   ============================================================================ */
"use strict";

/* ---------- DOM helpers ---------------------------------------------------- */
const $ = id => document.getElementById(id);
const qs = sel => document.querySelector(sel);
const qsa = sel => document.querySelectorAll(sel);

/* ---------- i18n ----------------------------------------------------------- *
 * Two-language dictionary (pt-BR / en).  All user-visible strings live here.
 * HTML elements are marked with `data-i18n` (textContent), `data-i18n-html`
 * (innerHTML, for elements that need <b>/<sub> markup) or
 * `data-i18n-attr-{title,placeholder}` (attributes).  Runtime strings —
 * toasts, dynamic warnings, button labels mid-operation — use `t(key)`.
 */
const I18N = {
  pt: {
    "nav.nanoeng": "NanoEng", "nav.research": "Pesquisa",
    "nav.tools": "Ferramentas", "nav.pubs": "Publicações",
    "header.theme":  "Alternar tema claro/escuro",

    "step1.title": "Etapa 1 — Estrutura",
    "step1.example": "Usar exemplo disponível",
    "step1.example.placeholder": "— selecionar —",
    "step1.load": "Carregar",
    "step1.or": "— ou envie seu próprio arquivo de estrutura —",
    "step1.drop.text": "<strong>Clique ou arraste</strong> um arquivo",
    "step1.drop.hint": "CIF · PDB · XYZ · POSCAR · XSF · LAMMPS · QE",
    "step1.info.lattice": "Rede", "step1.info.species": "Espécies",
    "step1.info.thetaSector": "setor θ",
    "step1.primitive": "🔬 Encontrar célula primitiva",
    "step1.primitive.finding": "Calculando…",
    "step1.cutoffs": "⚙ Cutoffs de ligação…",

    "step2.title": "Etapa 2 — Quiralidade",
    "step2.indices": "Índices quirais — clique no mapa ou digite abaixo",
    "step2.diameter": "Diâmetro", "step2.theta": "θ (quiral)",
    "step2.cellAtoms": "Átomos/cel.", "step2.strain": "Strain",

    "step3.title": "Etapa 3 — Parâmetros",
    "step3.reps": "Repetições axiais",
    "step3.vacuum": "Vácuo (Å)",
    "step3.rollin": "Roll inward",
    "step3.build": "Construir nanotubo",
    "step3.building": "Construindo…",
    "step3.undo": "↶ Desfazer", "step3.redo": "↷ Refazer",
    "step3.undo.title": "Desfazer última operação (Ctrl+Z)",
    "step3.redo.title": "Refazer (Ctrl+Y)",

    "step4.title": "Etapa 4 — Resultados &amp; ferramentas",
    "step4.chirality": "Quiralidade", "step4.atoms": "Átomos",
    "step4.diameter": "Diâmetro",    "step4.length": "Comprimento",
    "step4.strain": "Strain",        "step4.theta": "θ quiral",
    "step4.transforms": "Transformações avançadas",
    "step4.mwnt": "⧉ Multicamada…",
    "step4.deform": "↔ Deformação…",
    "step4.bundle": "⬡ Feixe…",
    "step4.analysis": "Análise &amp; relatórios",
    "step4.analysisBtn": "📈 Análise estrutural…",
    "step4.methodsBtn":  "📝 Texto de Methods…",
    "step4.dftBtn":      "⚛ Inputs DFT…",
    "step4.export": "Exportação",
    "step4.download": "⬇ Baixar estrutura…",
    "step4.batch":    "⊞ Lote (.zip)…",

    "map.title": "Mapa de quiralidade",
    "map.hint":  "clique para escolher (n,m)",
    "map.refresh": "Recalcular mapa",
    "map.dmin": "Diâmetro mín. (Å)",
    "map.dmax": "Diâmetro máx. (Å)",
    "map.strainMax": "Strain máx. (%)",
    "map.color": "Cor:",
    "map.atoms": "Átomos",
    "map.empty": "Carregue uma estrutura para gerar o mapa de quiralidade",
    "map.computing": "Calculando mapa de quiralidade…",
    "map.colorLabelAtoms":  "Átomos na célula",
    "map.colorLabelStrain": "Strain (%)",
    "map.legendClean":      "limpo",
    "map.legendSpurious":   "lig. espúrias",

    "viewer.title": "Visualizador 3D",
    "viewer.reset": "Reposicionar câmera",
    "viewer.style": "◉ Estilo", "viewer.style.title": "Alternar estilo",
    "viewer.bg": "◑ Fundo",    "viewer.bg.title": "Alternar fundo",
    "viewer.box": "□ Caixa",    "viewer.box.title": "Alternar caixa de simulação",
    "viewer.repz.title": "Replicar ao longo de z (apenas visual)",
    "viewer.repz.locked": "Desativado — torção quebra a periodicidade em Z",
    "viewer.empty": "Construa um nanotubo para visualizá-lo aqui",

    "footer.affil": "Universidade de Brasília — Laboratório de NanoEngenharia",
    "footer.cite":  "Se você usa o NTBuilder, por favor cite nosso paper.",
    "footer.github": "⬇ Versão local (GitHub)",

    "confirm.title": "Confirmar",
    "confirm.cancel": "Cancelar", "confirm.continue": "Continuar",
    "chain.title": "Operação encadeada — {op}",

    "batch.title": "⊞ Construção em lote",
    "batch.desc":  "Constrói vários nanotubos de uma vez e devolve um arquivo ZIP.",
    "batch.series": "Tipo de série",
    "batch.armchair": "Armchair (n, n)", "batch.zigzag": "Zigzag (n, 0)",
    "batch.all":      "Todos os pontos do mapa atual",
    "batch.nFrom": "n de", "batch.nTo": "n até",
    "batch.go":   "Construir &amp; baixar ZIP",
    "batch.going": "Construindo {n} tubos…",
    "batch.done":  "Lote concluído — {n} nanotubos.",
    "batch.fail":  "Falha no lote",

    "mwnt.title": "⧉ Nanotubo multicamada",
    "mwnt.desc":  "Constrói um nanotubo multi-parede escalonando o atual (n, m) por inteiros, de modo que todas as paredes compartilhem o mesmo período axial. Por design, esta operação reinicia a partir do (n, m) e ignora deformações já aplicadas.",
    "mwnt.walls": "Número de paredes",
    "mwnt.spacing": "Espaçamento interlamelar (Å)",
    "mwnt.rollin.hint": "inverte face interna/externa",
    "mwnt.go": "Construir MWNT",
    "mwnt.going": "Construindo…",

    "bundle.title": "⬡ Feixe de nanotubos",
    "bundle.desc":  "Replica o nanotubo atual em uma rede 2D formando uma supercélula periódica em bundle. Funciona sobre qualquer estrutura previamente construída (SWNT, MWNT ou deformada).",
    "bundle.geom": "Geometria",
    "bundle.linear": "Linear (2 tubos)",
    "bundle.triangle": "Triângulo (3 tubos)",
    "bundle.square4":  "Quadrado 2×2 (4 tubos)",
    "bundle.hex7":     "Hexagonal 1+6 (7 tubos)",
    "bundle.grid":     "Grade customizada (N × M)",
    "bundle.intertube": "Espaçamento entre tubos (Å)",
    "bundle.outerVac":  "Vácuo externo (Å)",
    "bundle.gridX": "Grade nx", "bundle.gridY": "Grade ny",
    "bundle.go": "Construir feixe", "bundle.going": "Construindo…",

    "deform.title": "↔ Deformação / Torção",
    "deform.desc":  "Aplica strain axial, strain radial e/ou torção ao nanotubo atual. A torção quebra a periodicidade axial: o resultado é um segmento finito em Z. A torção é aplicada sobre a supercélula que você está visualizando (controle ×z no cabeçalho do visualizador 3D). Para mudar o comprimento, feche este diálogo, ajuste o ×z e reabra.",
    "deform.strain":   "Strain axial (%)",
    "deform.twist":    "Taxa de torção (°/Å)",
    "deform.radial":   "Strain radial (%)",
    "deform.vacuum":   "Vácuo (Å)",
    "deform.zvac":     "Vácuo em Z (torção, Å)",
    "deform.go":   "Aplicar", "deform.going": "Aplicando…",

    "analysis.title": "📈 Análise estrutural",
    "analysis.desc":  "Histograma de comprimentos de ligação, caráter eletrônico por zone-folding (apenas grafeno) e classificação de simetria / grupo de linha do tubo.",
    "analysis.cutoff": "Cutoff de ligação (Å)",
    "analysis.run":    "Executar análise",
    "analysis.running": "Executando…",
    "analysis.bondStats": "Estatística de ligações",
    "analysis.bonds": "ligações",
    "analysis.mean": "média", "analysis.std": "σ",
    "analysis.min": "mín", "analysis.max": "máx",
    "analysis.pairs": "Pares de espécies:",
    "analysis.elec":  "Caráter eletrônico:",
    "analysis.sym":   "Simetria:",
    "analysis.xLabel": "Comprimento de ligação (Å)",
    "analysis.yLabel": "Contagem",

    "methods.title": "📝 Texto de Methods",
    "methods.desc":  "Parágrafo pronto para colar na seção Methods de um paper. A rede, quiralidade, caixa, vácuo, caráter eletrônico (quando aplicável) e qualquer deformação ativa são listados automaticamente.",
    "methods.cite": "BibTeX cite key",
    "methods.gen": "Gerar", "methods.copy": "Copiar",
    "methods.going": "Gerando…",
    "methods.copied": "Texto de Methods copiado.",

    "dft.title": "⚛ Arquivos de entrada DFT",
    "dft.desc":  "Gera arquivos de entrada com parâmetros sensatos para VASP, Quantum ESPRESSO, CP2K ou SIESTA a partir da estrutura atual. Revise e ajuste os parâmetros antes de lançar produção.",
    "dft.code": "Código",
    "dft.gen": "Gerar", "dft.copy": "Copiar",
    "dft.going": "Gerando…",
    "dft.copied": "Copiado: {name}",

    "cutoffs.title": "⚙ Cutoffs de ligação",
    "cutoffs.desc":  "Sobrescreve os cutoffs padrão por par de espécies usados na detecção de ligações espúrias induzidas por curvatura. Deixe vazio para usar os defaults internos.",
    "cutoffs.empty": "Sem sobrescritas — usando defaults internos.",
    "cutoffs.pair.placeholder": "Par (ex.: C-C, B-N)",
    "cutoffs.val.placeholder":  "cutoff (Å)",
    "cutoffs.add": "Adicionar",  "cutoffs.done": "Concluído",
    "cutoffs.invalid": "Informe o par (ex.: C-C) e o cutoff em Å.",

    "download.title": "⬇ Baixar estrutura",
    "download.desc":  "Escolha um formato. O arquivo é gerado pelo servidor a partir do estado atual da estrutura (incluindo deformações, multi-parede e feixe quando aplicáveis).",
    "download.close": "Fechar",
    "close": "Fechar",

    "rollin.help": "Inverte qual face do plano 2D fica voltada para dentro do tubo. Relevante apenas para redes heteronucleares (h-BN, MoS₂, MoSSe…).",

    // ─── Toasts / runtime strings ─────────────────────────────────────────
    "t.selectExample": "Selecione um exemplo primeiro.",
    "t.unsupportedFile": "Formato não suportado. Aceitos: CIF · PDB · XYZ · POSCAR/CONTCAR · XSF · LAMMPS · QE",
    "t.uploading": "Enviando…",
    "t.primAlreadyMsg": "A estrutura já está em forma primitiva.",
    "t.primFailed": "Falha em find_primitive_cell",
    "t.polarErr":   "Erro no mapa polar",
    "t.snapped":    "Simetria refinada",
    "t.undone": "Desfeito.", "t.redone": "Refeito.",
    "t.degen":  "(0, 0) é degenerado.",
    "t.buildOK.swnt": "SWNT construído",
    "t.buildOK.mwnt": "MWNT construído",
    "t.buildOK.bundle": "Feixe construído",
    "t.buildOK.deform": "Deformação construída",
    "t.buildFail":   "Falha ao construir",
    "t.loadStructFirst": "Carregue uma estrutura primeiro.",
    "t.buildFirst":      "Construa um nanotubo primeiro.",
    "t.opFail": "Falha em {label}",
    "t.opOK":   "{label}: {n} átomos",
    "t.analysisFail": "Falha na análise",
    "t.methodsFail":  "Falha ao gerar Methods",
    "t.dftFail":      "Falha ao gerar inputs DFT",
    "t.copyFail":     "Não foi possível copiar.",
    "t.deformProjExceeds": "Projeção ({n} átomos) excede o limite {max}. Reduza o ×z do visualizador antes de aplicar a torção.",
    "t.batchEmpty": "Nenhuma quiralidade válida.",
    "t.batchMax":   "Máximo de 100 quiralidades por lote.",
    "t.batchNeedMap": "Gere o mapa de quiralidade primeiro.",
    "t.batchOrder": "n inicial deve ser ≤ n final.",
    "t.deformApplied": "Aplicado",
    "t.deformChip":    "✦ Deformação aplicada",
    "t.repsLabel":     "Supercélula: {r} × {b} = {n} átomos",
    "t.wallsLabel":    "paredes",
    "t.wallPlanTitle": "Plano de paredes",
    "t.primFound": "Célula primitiva encontrada:",
    "t.primWas":   "eram",
    "t.primAtoms": "átomos",
    "t.primAlready": "Já é primitiva:",
    "t.projInfo":   "Projeção",
    "t.projOver":   "Excede o limite de {max} do /api/deform — reduza o ×z do visualizador.",

    // ─── Chain warnings ──────────────────────────────────────────────────
    "chain.mwnt.bundle":  "Você está prestes a construir um nanotubo multicamada cuja parede interna é um <b>feixe</b>. A rotina MWNT trata a estrutura atual como um único tubo — a geometria resultante provavelmente não é o que você espera.<br><br>Continuar mesmo assim?",
    "chain.mwnt.mwnt":    "A estrutura atual já é multicamada. Reexecutar Multicamada vai descartá-la e recomeçar a partir do (n, m) da parede interna.<br><br>Continuar?",
    "chain.mwnt.torsion": "Construir uma MWNT sobre um tubo torcido descarta a torção e recomeça a partir do (n, m) interno sem strain.<br><br>Continuar?",
    "chain.bundle.torsion": "Construir um feixe sobre um tubo torcido replica a geometria torcida, o que costuma quebrar a periodicidade do feixe.<br><br>Continuar?",
    "chain.axial.torsion":  "Aplicar strain axial sobre um tubo torcido produz densidade de torção não uniforme. Use apenas se souber o que está fazendo.<br><br>Continuar?",
    "chain.torsion.bundle": "Aplicar torção a um feixe quebra a periodicidade lateral. A geometria resultante é finita em Z.<br><br>Continuar?",
  },
  en: {
    "nav.nanoeng": "NanoEng", "nav.research": "Research",
    "nav.tools": "Tools", "nav.pubs": "Publications",
    "header.theme":  "Toggle light/dark theme",

    "step1.title": "Step 1 — Structure",
    "step1.example": "Use a bundled example",
    "step1.example.placeholder": "— select —",
    "step1.load": "Load",
    "step1.or": "— or upload your own structure file —",
    "step1.drop.text": "<strong>Click or drag</strong> a structure file",
    "step1.drop.hint": "CIF · PDB · XYZ · POSCAR · XSF · LAMMPS · QE",
    "step1.info.lattice": "Lattice", "step1.info.species": "Species",
    "step1.info.thetaSector": "θ sector",
    "step1.primitive": "🔬 Find primitive cell",
    "step1.primitive.finding": "Finding…",
    "step1.cutoffs": "⚙ Bond cutoffs…",

    "step2.title": "Step 2 — Chirality",
    "step2.indices": "Chiral indices — click the map or type below",
    "step2.diameter": "Diameter", "step2.theta": "θ (chiral)",
    "step2.cellAtoms": "Cell atoms", "step2.strain": "Strain",

    "step3.title": "Step 3 — Parameters",
    "step3.reps": "Axial repetitions",
    "step3.vacuum": "Vacuum (Å)",
    "step3.rollin": "Roll inward",
    "step3.build": "Build nanotube",
    "step3.building": "Building…",
    "step3.undo": "↶ Undo", "step3.redo": "↷ Redo",
    "step3.undo.title": "Undo last operation (Ctrl+Z)",
    "step3.redo.title": "Redo (Ctrl+Y)",

    "step4.title": "Step 4 — Results &amp; tools",
    "step4.chirality": "Chirality", "step4.atoms": "Atoms",
    "step4.diameter": "Diameter",   "step4.length": "Length",
    "step4.strain": "Strain",       "step4.theta": "θ chiral",
    "step4.transforms": "Advanced transformations",
    "step4.mwnt": "⧉ Multi-Wall…",
    "step4.deform": "↔ Deform / Torsion…",
    "step4.bundle": "⬡ Bundle…",
    "step4.analysis": "Analysis &amp; reports",
    "step4.analysisBtn": "📈 Structure analysis…",
    "step4.methodsBtn":  "📝 Methods text…",
    "step4.dftBtn":      "⚛ DFT inputs…",
    "step4.export": "Export",
    "step4.download": "⬇ Download structure…",
    "step4.batch":    "⊞ Batch (.zip)…",

    "map.title": "Chirality map",
    "map.hint":  "click to pick (n,m)",
    "map.refresh": "Recompute map",
    "map.dmin": "Min diameter (Å)",
    "map.dmax": "Max diameter (Å)",
    "map.strainMax": "Max strain (%)",
    "map.color": "Color:",
    "map.atoms": "Atoms",
    "map.empty": "Load a structure to generate the chirality map",
    "map.computing": "Computing chirality map…",
    "map.colorLabelAtoms":  "Atoms per unit cell",
    "map.colorLabelStrain": "Strain (%)",
    "map.legendClean":      "clean",
    "map.legendSpurious":   "spurious bonds",

    "viewer.title": "3D Viewer",
    "viewer.reset": "Reset camera",
    "viewer.style": "◉ Style", "viewer.style.title": "Cycle render style",
    "viewer.bg": "◑ BG",       "viewer.bg.title": "Toggle background",
    "viewer.box": "□ Box",     "viewer.box.title": "Toggle simulation box",
    "viewer.repz.title": "Replicate along z (display only)",
    "viewer.repz.locked": "Disabled — torsion breaks Z periodicity",
    "viewer.empty": "Build a nanotube to preview it here",

    "footer.affil": "University of Brasília — NanoEng Lab",
    "footer.cite":  "If you use NTBuilder, please cite our paper.",
    "footer.github": "⬇ Local install (GitHub)",

    "confirm.title": "Confirm",
    "confirm.cancel": "Cancel", "confirm.continue": "Continue",
    "chain.title": "Chained operation — {op}",

    "batch.title": "⊞ Batch build",
    "batch.desc":  "Build multiple nanotubes at once and download as a ZIP archive.",
    "batch.series": "Series type",
    "batch.armchair": "Armchair (n, n)", "batch.zigzag": "Zigzag (n, 0)",
    "batch.all":      "All points on current map",
    "batch.nFrom": "n from", "batch.nTo": "n to",
    "batch.go":   "Build &amp; download ZIP",
    "batch.going": "Building {n} tubes…",
    "batch.done":  "Batch complete — {n} nanotubes.",
    "batch.fail":  "Batch failed",

    "mwnt.title": "⧉ Multi-Wall Nanotube",
    "mwnt.desc":  "Build a multi-walled nanotube by integer scaling of the current (n, m) so every wall shares the same axial period. By design this routine restarts from (n, m) and ignores any deformation chain.",
    "mwnt.walls": "Number of walls",
    "mwnt.spacing": "Interlayer spacing (Å)",
    "mwnt.rollin.hint": "invert inside/outside face",
    "mwnt.go": "Build MWNT",
    "mwnt.going": "Building…",

    "bundle.title": "⬡ Nanotube bundle",
    "bundle.desc":  "Replicate the current nanotube on a 2D lattice to form a periodic bundle supercell. Works on top of any previously built structure (SWNT, MWNT or deformed tube).",
    "bundle.geom": "Geometry",
    "bundle.linear": "Linear (2 tubes)",
    "bundle.triangle": "Triangle (3 tubes)",
    "bundle.square4":  "Square 2×2 (4 tubes)",
    "bundle.hex7":     "Hexagonal 1+6 (7 tubes)",
    "bundle.grid":     "Custom grid (N × M)",
    "bundle.intertube": "Inter-tube spacing (Å)",
    "bundle.outerVac":  "Outer vacuum (Å)",
    "bundle.gridX": "Grid nx", "bundle.gridY": "Grid ny",
    "bundle.go": "Build bundle", "bundle.going": "Building…",

    "deform.title": "↔ Deform / Torsion",
    "deform.desc":  "Apply axial strain, radial strain and / or torsion to the current nanotube. Torsion breaks axial periodicity: the result is a finite Z segment. The twist is applied to the supercell currently being shown (×z control in the 3D viewer header). To change the length, close this dialog, adjust ×z and reopen.",
    "deform.strain":   "Axial strain (%)",
    "deform.twist":    "Torsion rate (°/Å)",
    "deform.radial":   "Radial strain (%)",
    "deform.vacuum":   "Vacuum (Å)",
    "deform.zvac":     "Z vacuum (torsion, Å)",
    "deform.go":   "Apply", "deform.going": "Applying…",

    "analysis.title": "📈 Structure analysis",
    "analysis.desc":  "Bond-length histogram, zone-folding electronic character (graphene only) and tube symmetry / line-group classification.",
    "analysis.cutoff": "Bond cutoff (Å)",
    "analysis.run":    "Run analysis",
    "analysis.running": "Running…",
    "analysis.bondStats": "Bond statistics",
    "analysis.bonds": "bonds",
    "analysis.mean": "mean", "analysis.std": "σ",
    "analysis.min": "min", "analysis.max": "max",
    "analysis.pairs": "Species pairs:",
    "analysis.elec":  "Electronic character:",
    "analysis.sym":   "Symmetry:",
    "analysis.xLabel": "Bond length (Å)",
    "analysis.yLabel": "Count",

    "methods.title": "📝 Methods text",
    "methods.desc":  "Ready-to-paste paragraph for the Methods section of a paper. The lattice, chirality, box, vacuum, electronic character (when applicable) and any active deformation are listed automatically.",
    "methods.cite": "BibTeX cite key",
    "methods.gen": "Generate", "methods.copy": "Copy",
    "methods.going": "Generating…",
    "methods.copied": "Methods text copied.",

    "dft.title": "⚛ DFT input files",
    "dft.desc":  "Generate sensible starting input files for VASP, Quantum ESPRESSO, CP2K or SIESTA from the current structure. Review the parameters before launching production.",
    "dft.code": "Code",
    "dft.gen": "Generate", "dft.copy": "Copy",
    "dft.going": "Generating…",
    "dft.copied": "Copied: {name}",

    "cutoffs.title": "⚙ Bond cutoffs",
    "cutoffs.desc":  "Override the default per-pair bond cutoffs used by the curvature-induced spurious-bond check. Leave empty to use built-in defaults.",
    "cutoffs.empty": "No overrides — using built-in defaults.",
    "cutoffs.pair.placeholder": "Pair (e.g. C-C, B-N)",
    "cutoffs.val.placeholder":  "cutoff (Å)",
    "cutoffs.add": "Add", "cutoffs.done": "Done",
    "cutoffs.invalid": "Provide a pair (e.g. C-C) and a cutoff in Å.",

    "download.title": "⬇ Download structure",
    "download.desc":  "Pick a format. The file is generated server-side from the current structure (including any active deformation, multi-wall and bundle).",
    "download.close": "Close",
    "close": "Close",

    "rollin.help": "Inverts which face of the 2D plane points inward in the tube. Relevant only for heteronuclear sheets (h-BN, MoS₂, MoSSe…).",

    "t.selectExample": "Select an example first.",
    "t.unsupportedFile": "Unsupported file. Accepted: CIF · PDB · XYZ · POSCAR/CONTCAR · XSF · LAMMPS · QE",
    "t.uploading": "Uploading…",
    "t.primAlreadyMsg": "Structure is already primitive.",
    "t.primFailed": "find_primitive_cell failed",
    "t.polarErr":   "Polar map error",
    "t.snapped":    "Symmetry snapped",
    "t.undone": "Undone.", "t.redone": "Redone.",
    "t.degen":  "(0, 0) is degenerate.",
    "t.buildOK.swnt": "SWNT built",
    "t.buildOK.mwnt": "MWNT built",
    "t.buildOK.bundle": "Bundle built",
    "t.buildOK.deform": "Deformation applied",
    "t.buildFail":   "Build failed",
    "t.loadStructFirst": "Load a structure first.",
    "t.buildFirst":      "Build a nanotube first.",
    "t.opFail": "{label} failed",
    "t.opOK":   "{label}: {n} atoms",
    "t.analysisFail": "Analysis failed",
    "t.methodsFail":  "Methods generation failed",
    "t.dftFail":      "DFT inputs failed",
    "t.copyFail":     "Could not copy to clipboard.",
    "t.deformProjExceeds": "Projection ({n} atoms) exceeds the {max} limit. Reduce the viewer's ×z before applying the torsion.",
    "t.batchEmpty": "No valid chiralities.",
    "t.batchMax":   "Maximum 100 chiralities per batch.",
    "t.batchNeedMap": "Generate the chirality map first.",
    "t.batchOrder": "n_from must be ≤ n_to.",
    "t.deformApplied": "Applied",
    "t.deformChip":    "✦ Deformation applied",
    "t.repsLabel":     "Supercell: {r} × {b} = {n} atoms",
    "t.wallsLabel":    "walls",
    "t.wallPlanTitle": "Wall plan",
    "t.primFound": "Primitive cell found:",
    "t.primWas":   "was",
    "t.primAtoms": "atoms",
    "t.primAlready": "Already primitive:",
    "t.projInfo":   "Projection",
    "t.projOver":   "Exceeds the {max} limit of /api/deform — reduce the viewer's ×z.",

    "chain.mwnt.bundle":  "You are about to build a multi-walled nanotube whose inner wall is a <b>bundle</b>. The MWNT pipeline treats the current structure as a single tube — the geometry probably isn't what you expect.<br><br>Continue anyway?",
    "chain.mwnt.mwnt":    "The current structure is already multi-walled. Re-running Multi-Wall will discard it and start from the inner (n, m) again.<br><br>Continue?",
    "chain.mwnt.torsion": "Building an MWNT on top of a torqued tube discards the torsion and restarts from the unstrained inner (n, m).<br><br>Continue?",
    "chain.bundle.torsion": "Building a bundle on top of a torqued tube replicates the twisted geometry, which typically breaks the bundle's periodicity.<br><br>Continue anyway?",
    "chain.axial.torsion":  "Axial strain on top of a torqued tube produces a non-uniform twist density. Use only if you know what you're doing.<br><br>Continue?",
    "chain.torsion.bundle": "Applying torsion to a bundle breaks its lateral periodicity. The resulting geometry is finite along Z.<br><br>Continue?",
  },
};

let _LANG = "pt";
function t(key, vars) {
  const dict = I18N[_LANG] || I18N.pt;
  let s = dict[key];
  if (s === undefined) s = (I18N.pt[key] !== undefined) ? I18N.pt[key] : key;
  if (vars) {
    for (const [k, v] of Object.entries(vars))
      s = s.replace(new RegExp(`\\{${k}\\}`, "g"), String(v));
  }
  return s;
}

function applyLang(lang) {
  if (!I18N[lang]) lang = "pt";
  _LANG = lang;
  document.documentElement.lang = lang === "en" ? "en" : "pt-BR";
  try { localStorage.setItem("ntbuilder-lang", lang); } catch (_) {}
  qsa("[data-i18n]").forEach(el => {
    const k = el.dataset.i18n;
    const v = I18N[lang][k];
    if (v !== undefined) el.textContent = v;
  });
  qsa("[data-i18n-html]").forEach(el => {
    const k = el.dataset.i18nHtml;
    const v = I18N[lang][k];
    if (v !== undefined) el.innerHTML = v;
  });
  qsa("[data-i18n-title]").forEach(el => {
    const k = el.dataset.i18nTitle;
    const v = I18N[lang][k];
    if (v !== undefined) el.title = v;
  });
  qsa("[data-i18n-placeholder]").forEach(el => {
    const k = el.dataset.i18nPlaceholder;
    const v = I18N[lang][k];
    if (v !== undefined) el.placeholder = v;
  });
  // Sync the segmented toggle visual state
  qsa(".lang-toggle button").forEach(b => {
    const on = b.dataset.lang === lang;
    b.classList.toggle("active", on);
    b.setAttribute("aria-pressed", on ? "true" : "false");
  });
  // Refresh dynamic UI that depends on language
  if (state && state.summary) updateAfterBuild(false);
  if (typeof updateDeformProjection === "function") updateDeformProjection();
  if (typeof renderCutoffs === "function") renderCutoffs();
  if (state && state.polarData) renderPolar(state.polarData);
}

document.addEventListener("click", e => {
  const b = e.target.closest(".lang-toggle button[data-lang]");
  if (b) applyLang(b.dataset.lang);
});

/* ---------- Global state --------------------------------------------------- */
const state = {
  // structure / examples
  fileId:    null,
  example:   null,
  polarData: null,
  bondCutoffs: {},    // per-pair Å overrides for /api/build

  // current job
  jobId:    null,
  built:    false,
  xyzStr:   null,
  box:      null,
  kind:     null,     // "swnt" | "mwnt" | "bundle" | "axial" | "radial" | "torsion"
  summary:  null,     // {n, m, diameter, length, n_atoms, strain, theta_deg}
  walls:    null,     // wall plan from /api/mwnt
  warning:  null,
  torsionApplied: false,
  repsAtBuild: 1,     // axial-reps used in the latest deform that may bake torsion
  atomsBase: null,    // atoms in 1 axial period (used to display reps × atoms_base)
  deformDesc: null,

  // history
  undoStack: [],
  redoStack: [],

  // ui
  colorMode: "strain",
  viewer:    null,
  viewStyle: 0,
  viewBg:    null,    // override BG (null = follow theme)
  showBox:   true,
};

/* ---------- API helpers ---------------------------------------------------- */
async function api(method, path, body = null, isFormData = false) {
  const opts = { method };
  if (body) {
    if (isFormData) {
      opts.body = body;
    } else {
      opts.headers = { "Content-Type": "application/json" };
      opts.body = JSON.stringify(body);
    }
  }
  const res = await fetch(path, opts);
  if (!res.ok) {
    let detail = `HTTP ${res.status}`;
    try {
      const data = await res.json();
      detail = data.detail || JSON.stringify(data);
    } catch (_) { /* ignore */ }
    throw new Error(detail);
  }
  return res;
}
const apiJSON = async (...args) => (await api(...args)).json();

/* ---------- Toast ---------------------------------------------------------- */
function toast(msg, type = "info", ms = 3500) {
  const t = document.createElement("div");
  t.className = `toast ${type}`;
  t.textContent = msg;
  $("toast-container").appendChild(t);
  setTimeout(() => t.remove(), ms);
}

/* ---------- Theme ---------------------------------------------------------- */
function applyTheme(theme) {
  document.documentElement.setAttribute("data-theme", theme);
  $("theme-toggle-btn").textContent = theme === "dark" ? "☀" : "🌙";
  try { localStorage.setItem("ntbuilder-theme", theme); } catch (_) {}
  if (state.viewer) {
    state.viewer.setBackgroundColor(getViewerBg());
    state.viewer.render();
  }
  if (state.polarData) renderPolar(state.polarData);
}
$("theme-toggle-btn").addEventListener("click", () => {
  const cur = document.documentElement.getAttribute("data-theme") || "light";
  applyTheme(cur === "dark" ? "light" : "dark");
});
(function () {
  let saved = "light";
  try { saved = localStorage.getItem("ntbuilder-theme") || "light"; } catch (_) {}
  applyTheme(saved);
})();

function getViewerBg() {
  if (state.viewBg) return state.viewBg;
  const dark = document.documentElement.getAttribute("data-theme") === "dark";
  return dark ? "#060B17" : "#F8FAFC";
}

function plotColors() {
  const dark = document.documentElement.getAttribute("data-theme") === "dark";
  return {
    grid:   dark ? "#1F2A44" : "#E2E8F0",
    text:   dark ? "#94A3B8" : "#64748B",
    radial: dark ? "rgba(255,255,255,.18)" : "rgba(15,23,42,.20)",
  };
}

/* ---------- Examples ------------------------------------------------------- */
async function loadExamples() {
  try {
    const list = await apiJSON("GET", "api/examples");
    list.forEach(e => {
      const opt = document.createElement("option");
      opt.value = e.filename;
      opt.textContent = e.name;
      $("example-select").appendChild(opt);
    });
  } catch (err) {
    console.warn("examples: ", err);
  }
}
$("example-load-btn").addEventListener("click", () => {
  const val = $("example-select").value;
  if (!val) { toast(t("t.selectExample"), "err"); return; }
  state.fileId = null;
  state.example = val;
  setCifStatus("ok", val.replace(/\.[^.]+$/, "").replace(/_/g, " "));
  $("btn-prim").disabled = false;
  $("prim-result").classList.add("hidden");
  runPolar();
});

/* ---------- File upload ---------------------------------------------------- */
const _ACCEPTED_EXTS = new Set([
  ".cif", ".pdb", ".xyz",
  ".poscar", ".contcar", ".vasp",
  ".xsf", ".lammps", ".data",
  ".in", ".pwi",
]);
const _ACCEPTED_NAMES = new Set(["poscar", "contcar"]);

const dropZone = $("drop-zone");
dropZone.addEventListener("dragover", e => {
  e.preventDefault(); dropZone.classList.add("dragging");
});
dropZone.addEventListener("dragleave", () => dropZone.classList.remove("dragging"));
dropZone.addEventListener("drop", e => {
  e.preventDefault();
  dropZone.classList.remove("dragging");
  if (e.dataTransfer.files[0]) handleFile(e.dataTransfer.files[0]);
});
$("cif-file-input").addEventListener("change", e => {
  if (e.target.files[0]) handleFile(e.target.files[0]);
});

async function handleFile(file) {
  const name = file.name.toLowerCase();
  const ext = "." + name.split(".").pop();
  const stem = name.split(".")[0];
  if (!_ACCEPTED_EXTS.has(ext) && !_ACCEPTED_NAMES.has(stem)) {
    toast(t("t.unsupportedFile"), "err", 5000);
    return;
  }
  setCifStatus("info", t("t.uploading"));
  try {
    const fd = new FormData();
    fd.append("file", file);
    const res = await apiJSON("POST", "api/upload", fd, true);
    state.fileId = res.file_id;
    state.example = null;
    setCifStatus("ok", file.name);
    $("btn-prim").disabled = false;
    $("prim-result").classList.add("hidden");
    runPolar();
  } catch (err) {
    setCifStatus("err", err.message);
  }
}

function setCifStatus(kind, msg) {
  const el = $("cif-status");
  el.className = "";
  if (kind === "info") {
    el.classList.add("info");
    el.innerHTML = `<span class="spinner" style="width:12px;height:12px;"></span> ${msg}`;
  } else if (kind === "ok") {
    el.classList.add("ok");
    el.innerHTML = `✓ ${msg}`;
    $("build-btn").disabled = false;
    $("batch-open-btn").disabled = false;
  } else {
    el.classList.add("err");
    el.innerHTML = `✗ ${msg}`;
  }
}

/* ---------- Primitive cell ------------------------------------------------- */
$("btn-prim").addEventListener("click", async () => {
  if (!state.fileId && !state.example) return;
  const btn = $("btn-prim");
  btn.disabled = true;
  const old = btn.innerHTML;
  btn.innerHTML = `<span class="spinner" style="width:12px;height:12px;margin-right:6px;"></span>${t("step1.primitive.finding")}`;
  try {
    const qry = state.fileId
      ? `file_id=${encodeURIComponent(state.fileId)}`
      : `example=${encodeURIComponent(state.example)}`;
    const res = await apiJSON("POST", `api/primitive?${qry}`);
    const changed = res.n_atoms_prim !== res.n_atoms_orig
                 || Math.abs(res.a_prim - res.a_orig) > 0.001;
    const out = $("prim-result");
    if (changed) {
      if (res.file_id) { state.fileId = res.file_id; state.example = null; }
      out.innerHTML =
        `<b>${t("t.primFound")}</b> ${res.n_atoms_prim} ${t("t.primAtoms")} `
        + `(${t("t.primWas")} ${res.n_atoms_orig}) · a=${res.a_prim} Å · ${res.lattice_type}<br>`
        + `<span style="color:var(--c-muted);font-size:.7rem;">${res.description}</span>`;
      toast(`${t("t.primFound").replace(":", "")} ${res.n_atoms_prim} ${t("t.primAtoms")} (${res.lattice_type})`, "ok");
      runPolar();
    } else {
      out.innerHTML = `<b>${t("t.primAlready")}</b> ${res.n_atoms_orig} ${t("t.primAtoms")} · a=${res.a_orig} Å`;
      toast(t("t.primAlreadyMsg"), "info");
    }
    out.classList.remove("hidden");
  } catch (err) {
    toast(t("t.primFailed") + ": " + err.message, "err", 5000);
  } finally {
    btn.textContent = "";
    btn.innerHTML = t("step1.primitive");   // restaura idioma corrente
    btn.disabled = false;
  }
});

/* ---------- Polar map ------------------------------------------------------ */
async function runPolar() {
  const empty = $("polar-empty");
  empty.style.display = "flex";
  empty.innerHTML = `<div class="spinner"></div><span>${t("map.computing")}</span>`;

  const dMax = parseFloat($("inp-dmax").value) || 25;
  const useStrainFilter = $("strain-filter-chk").checked;
  // ``roll_inward`` from the sidebar feeds the spurious-bond check on
  // the server side: flipping the toggle swaps the concave / convex face
  // of a buckled monolayer, so chiralities that were clean may become
  // marked × and vice-versa.  The sidebar toggle is the canonical source
  // of truth even for the polar map (so that the X markers always match
  // the build that would actually be produced if the user clicked a
  // point right now).
  const rollInChk = $("inp-rollin");
  const rollInward = rollInChk ? rollInChk.checked : false;
  const body = {
    max_diameter: dMax,
    n_max: 60,
    roll_inward: rollInward,
    ...(useStrainFilter ? { strain_max: parseFloat($("inp-strain-max").value) || 5 } : {}),
    ...(state.fileId  ? { file_id: state.fileId  } : {}),
    ...(state.example ? { example: state.example } : {}),
  };

  try {
    const data = await apiJSON("POST", "api/polar", body);
    state.polarData = data;
    empty.style.display = "none";
    renderPolar(data);
    $("polar-info").textContent = `${data.points.length} pts · a=${data.a} Å · ${data.lattice_type}`;
    updateStructInfo(data);
    if (data.snap_desc) toast(`${t("t.snapped")}: ${data.snap_desc}`, "info", 5000);
  } catch (err) {
    empty.innerHTML = `<div class="big-icon">⚠</div><span>${err.message}</span>`;
    toast(t("t.polarErr") + ": " + err.message, "err");
  }
}

function updateStructInfo(data) {
  $("si-lattice").textContent = data.lattice_type || "—";
  $("si-ab").textContent      = data.a === data.b ? `${data.a} Å` : `${data.a} / ${data.b} Å`;
  $("si-gamma").textContent   = `${data.gamma_deg}°`;
  $("si-species").textContent = (data.species || []).join(", ") || "—";
  $("si-dmin").textContent    = `${data.d_min} Å`;
  $("si-theta").textContent   = `${data.theta_max.toFixed(1)}°`;
  $("struct-info").classList.remove("hidden");

  // Roll inward só faz diferença quando há mais de uma espécie distinguível.
  // Para sistemas monatômicos (graphene, irida-graphene, biphenylene...) o
  // toggle é redundante e atrapalha a UI — escondemos.
  const nSpecies = new Set(data.species || []).size;
  const showRollIn = nSpecies > 1;
  ["rollin-row", "mwnt-rollin-cell", "batch-rollin-row"].forEach(id => {
    const el = $(id);
    if (!el) return;
    el.style.display = showRollIn ? "" : "none";
    if (!showRollIn) {
      // Resetar valor para evitar enviar true acidentalmente em payloads futuros
      const cb = el.querySelector('input[type="checkbox"]');
      if (cb) cb.checked = false;
    }
  });
}

$("polar-refresh-btn").addEventListener("click", () => {
  if (state.fileId || state.example) runPolar();
});

$("strain-filter-chk").addEventListener("change", () => {
  $("inp-strain-max").disabled = !$("strain-filter-chk").checked;
  if (state.fileId || state.example) runPolar();
});
$("inp-strain-max").addEventListener("change", () => {
  if ($("strain-filter-chk").checked && (state.fileId || state.example)) runPolar();
});

// "Roll inward" toggle (Step 3 sidebar) feeds the spurious-bond check on
// /api/polar: flipping it swaps which face of a Janus monolayer faces the
// concave side, so the set of × markers on the map changes too.  We
// re-fetch instead of merely re-rendering because the actual list of
// spurious pairs has to be recomputed on the server.
const _rollinSidebar = $("inp-rollin");
if (_rollinSidebar) {
  _rollinSidebar.addEventListener("change", () => {
    if (state.fileId || state.example) runPolar();
  });
}

["change", "input"].forEach(ev => {
  $("inp-dmax").addEventListener(ev, () => {
    if (state.fileId || state.example) {
      clearTimeout(window._dmaxTimer);
      window._dmaxTimer = setTimeout(runPolar, 300);
    }
  });
  $("inp-dmin").addEventListener(ev, () => {
    if (state.polarData) renderPolar(state.polarData);
  });
});

/* Color segmented switch */
qsa("#color-seg button").forEach(b => {
  b.addEventListener("click", () => {
    qsa("#color-seg button").forEach(x => x.classList.remove("active"));
    b.classList.add("active");
    state.colorMode = b.dataset.color;
    if (state.polarData) renderPolar(state.polarData);
  });
});

function getFilteredPoints(data) {
  const dMin = parseFloat($("inp-dmin").value) || 0;
  const dMax = parseFloat($("inp-dmax").value) || 9999;
  return data.points.filter(p => p.diameter >= dMin && p.diameter <= dMax);
}

function renderPolar(data) {
  const pts = getFilteredPoints(data);
  if (!pts.length) return;

  const c = plotColors();
  const mode = state.colorMode;
  const colorVals = pts.map(p => mode === "n_atoms" ? p.n_atoms : p.strain);
  const colorLabel = mode === "n_atoms" ? t("map.colorLabelAtoms") : t("map.colorLabelStrain");
  const cMin = Math.min(...colorVals), cMax = Math.max(...colorVals) || 1e-6;
  const dmax = data.dmax, thetaMax = data.theta_max;
  const tmRad = thetaMax * Math.PI / 180;

  function arc(r, nPts = 80) {
    const ax = [], ay = [];
    for (let i = 0; i <= nPts; i++) {
      const t = (i / nPts) * tmRad;
      ax.push(r * Math.cos(t)); ay.push(r * Math.sin(t));
    }
    return { ax, ay };
  }
  function concat(...segs) {
    const x = [], y = [];
    segs.forEach((s, i) => {
      if (i > 0) { x.push(null); y.push(null); }
      x.push(...s.x); y.push(...s.y);
    });
    return { x, y };
  }

  const dStep = dmax <= 15 ? 5 : dmax <= 30 ? 5 : 10;
  const arcSegs = [], dLabels = [];
  for (let d = dStep; d < dmax; d += dStep) {
    const { ax, ay } = arc(d);
    arcSegs.push({ x: ax, y: ay });
    dLabels.push({ x: d, y: -dmax * 0.03, text: `${d}`, showarrow: false,
                   font: { size: 8, color: c.text }, xanchor: "center", yanchor: "top" });
  }
  const gridArcs = concat(...arcSegs);
  const gridArcTrace = {
    type: "scatter", mode: "lines",
    x: gridArcs.x, y: gridArcs.y,
    line: { color: c.grid, width: 0.7 },
    hoverinfo: "skip", showlegend: false,
  };

  const tStep = thetaMax <= 30 ? 10 : thetaMax <= 60 ? 15 : 30;
  const angSegs = [], tLabels = [];
  for (let th = tStep; th < thetaMax; th += tStep) {
    const t = th * Math.PI / 180;
    angSegs.push({ x: [0, dmax * Math.cos(t)], y: [0, dmax * Math.sin(t)] });
    tLabels.push({ x: (dmax + 1.2) * Math.cos(t), y: (dmax + 1.2) * Math.sin(t),
                   text: `${th}°`, showarrow: false,
                   font: { size: 8, color: c.text }, xanchor: "center", yanchor: "middle" });
  }
  const gridAngTrace = angSegs.length ? {
    type: "scatter", mode: "lines",
    ...concat(...angSegs),
    line: { color: c.grid, width: 0.7 },
    hoverinfo: "skip", showlegend: false,
  } : { type: "scatter", mode: "lines", x: [], y: [], hoverinfo: "skip", showlegend: false };

  const { ax: bax, ay: bay } = arc(dmax);
  const boundaryTrace = {
    type: "scatter", mode: "lines",
    x: bax, y: bay,
    line: { color: c.radial, width: 1.5 },
    hoverinfo: "skip", showlegend: false,
  };
  const boundaryLinesTrace = {
    type: "scatter", mode: "lines",
    ...concat(
      { x: [0, dmax * 1.02], y: [0, 0] },
      { x: [0, dmax * 1.02 * Math.cos(tmRad)], y: [0, dmax * 1.02 * Math.sin(tmRad)] },
    ),
    line: { color: c.radial, width: 1, dash: "dot" },
    hoverinfo: "skip", showlegend: false,
  };

  // Curvature-induced spurious bonds (e.g. Se-Se in MoSSe small-D tubes)
  // are surfaced as an "x" marker on the affected chiralities.  This is
  // informational only: the user may still click and build the tube.
  // Plotly accepts ``marker.symbol`` / ``marker.size`` as per-point
  // arrays, so we keep a single trace (preserves the index → 5 of
  // ``selectedTrace`` that downstream Plotly.restyle calls depend on).
  const markerSymbols = pts.map(p =>
    (p.spurious && p.spurious.length) ? "x" : "circle");
  const markerSizes = pts.map(p =>
    (p.spurious && p.spurious.length) ? 10 : 7);

  const trace = {
    type: "scatter", mode: "markers",
    x: pts.map(p => p.x), y: pts.map(p => p.y),
    marker: {
      symbol: markerSymbols,
      size:   markerSizes,
      color: colorVals,
      colorscale: "Viridis",
      cmin: cMin, cmax: cMax,
      showscale: true,
      colorbar: {
        orientation: "h", y: -0.14, x: 0.5, xanchor: "center", yanchor: "top",
        len: 0.75, thickness: 11,
        title: { text: colorLabel, side: "bottom", font: { size: 9, color: c.text } },
        tickfont: { size: 9, color: c.text },
        bgcolor: "transparent", bordercolor: "transparent", outlinecolor: "transparent",
      },
      line: { width: 0.5, color: "rgba(120,130,150,.25)" },
    },
    text: pts.map(p => {
      const sp = (p.spurious && p.spurious.length)
        ? `<br><b>spurious bonds:</b> ${p.spurious.join(", ")}` : "";
      return `(${p.n},${p.m})<br>D=${p.diameter} Å<br>θ=${p.theta_deg.toFixed(1)}°<br>`
           + `atoms=${p.n_atoms}<br>strain=${p.strain.toFixed(4)}%${sp}`;
    }),
    hovertemplate: "%{text}<extra></extra>",
    customdata: pts,
  };

  const selectedTrace = {
    type: "scatter", mode: "markers",
    x: [], y: [],
    marker: { size: 14, color: "rgba(0,0,0,0)", line: { width: 2.5, color: "#F59E0B" } },
    hoverinfo: "skip", showlegend: false, name: "selected",
  };

  const nStep = Math.max(5, Math.round(dmax / 5 / 5) * 5);
  const zigzagLabels = data.points
    .filter(p => p.m === 0 && p.n > 0 && p.n % nStep === 0)
    .map(p => ({ x: p.x, y: -dmax * 0.055, text: `n=${p.n}`, showarrow: false,
                 font: { size: 8, color: c.text }, xanchor: "center", yanchor: "top" }));
  const armLabels = data.points
    .filter(p => p.n === p.m && p.n > 0 && p.n % nStep === 0)
    .map(p => {
      const off = dmax * 0.06;
      return { x: p.x - off * Math.sin(tmRad), y: p.y + off * Math.cos(tmRad),
               text: `n=${p.n}`, showarrow: false,
               font: { size: 8, color: c.text }, xanchor: "right", yanchor: "middle" };
    });

  // A "legend" describing what each marker means.  We only render it for
  // structures that can develop curvature-induced spurious bonds (i.e.
  // those with buckling) — for flat lattices both symbols would have the
  // same meaning and the legend would be noise.  Whether the legend is
  // shown is driven by the *presence of the spurious field* on the
  // points (set server-side only when the structure has buckling),
  // independent of whether any individual chirality is currently marked.
  const anyHasSpuriousField = pts.some(p => Array.isArray(p.spurious));
  const xLegend = dmax * 1.04;
  const yLegendTop = dmax * 0.92 * Math.sin(tmRad);
  const yLegendStep = dmax * 0.06;
  const spuriousLegend = anyHasSpuriousField ? [
    {
      x: xLegend, y: yLegendTop,
      text: `● <span style="font-size:8px">${t("map.legendClean")}</span>`,
      showarrow: false,
      font: { size: 11, color: c.text },
      xanchor: "left", yanchor: "middle",
    },
    {
      x: xLegend, y: yLegendTop - yLegendStep,
      text: `✕ <span style="font-size:8px">${t("map.legendSpurious")}</span>`,
      showarrow: false,
      font: { size: 11, color: c.text },
      xanchor: "left", yanchor: "middle",
    },
  ] : [];

  const annotations = [
    { x: dmax * 1.05, y: -dmax * 0.03, text: "D (Å)", showarrow: false,
      font: { size: 8, color: c.text }, xanchor: "left", yanchor: "top" },
    { x: dmax * 1.04, y: 0, text: "<b>Zigzag</b> (m=0)", showarrow: false,
      font: { size: 9, color: c.text }, xanchor: "left", yanchor: "middle" },
    { x: dmax * 1.04 * Math.cos(tmRad), y: dmax * 1.04 * Math.sin(tmRad),
      text: "<b>Armchair</b> (n=m)", showarrow: false,
      font: { size: 9, color: c.text }, xanchor: "left", yanchor: "bottom" },
    { x: 0, y: -dmax * 0.03, text: "0", showarrow: false,
      font: { size: 8, color: c.text }, xanchor: "center", yanchor: "top" },
    ...dLabels, ...tLabels, ...zigzagLabels, ...armLabels, ...spuriousLegend,
  ];

  const yMax = dmax * Math.sin(tmRad);
  const layout = {
    paper_bgcolor: "transparent", plot_bgcolor: "transparent",
    margin: { l: 20, r: 20, t: 8, b: 60 },
    xaxis: { range: [-dmax * 0.05, dmax * 1.18], color: c.text,
             gridcolor: "transparent", zerolinecolor: "transparent",
             tickfont: { size: 8, color: c.text }, showticklabels: false },
    yaxis: { range: [-dmax * 0.08, yMax * 1.15], color: c.text,
             gridcolor: "transparent", zerolinecolor: "transparent",
             tickfont: { size: 8, color: c.text }, showticklabels: false,
             scaleanchor: "x", scaleratio: 1 },
    annotations,
    font: { family: "Inter, system-ui, sans-serif", size: 9, color: c.text },
    showlegend: false,
  };

  const div = $("polar-plot");
  Plotly.react(div,
    [gridArcTrace, gridAngTrace, boundaryTrace, boundaryLinesTrace, trace, selectedTrace],
    layout, { responsive: true, displayModeBar: false });
  div.on("plotly_click", evt => {
    const pt = evt.points.find(p => p.data.customdata);
    if (!pt) return;
    selectChirality(pt.data.customdata[pt.pointIndex]);
  });
}

function selectChirality(d) {
  $("inp-n").value = d.n; $("inp-m").value = d.m;
  $("cp-d").textContent      = `${d.diameter} Å`;
  $("cp-theta").textContent  = `${d.theta_deg.toFixed(2)}°`;
  $("cp-atoms").textContent  = d.n_atoms;
  $("cp-strain").textContent = `${d.strain.toFixed(4)} %`;
  $("chiral-props").classList.remove("hidden");
  Plotly.restyle("polar-plot", { x: [[d.x]], y: [[d.y]] }, [5]);
}

function onNMChange() {
  if (!state.polarData) return;
  const n = parseInt($("inp-n").value, 10), m = parseInt($("inp-m").value, 10);
  const hit = state.polarData.points.find(p => p.n === n && p.m === m);
  if (hit) selectChirality(hit);
  else {
    Plotly.restyle("polar-plot", { x: [[]], y: [[]] }, [5]);
    $("chiral-props").classList.add("hidden");
  }
}
$("inp-n").addEventListener("input", onNMChange);
$("inp-m").addEventListener("input", onNMChange);

/* ---------- Reps display & torsion-hides-reps ------------------------------ */
function updateRepsInfo() {
  const reps = parseInt($("inp-repeat").value, 10) || 1;
  const el = $("reps-info");
  if (state.atomsBase && reps > 1) {
    el.textContent = t("t.repsLabel", {
      r: reps, b: state.atomsBase,
      n: (reps * state.atomsBase).toLocaleString(),
    });
    el.classList.remove("hidden");
  } else {
    el.classList.add("hidden");
  }
}
$("inp-repeat").addEventListener("input", updateRepsInfo);

function setRepsVisible(visible) {
  const row = $("reps-row");
  if (!row) return;
  if (visible) {
    row.style.display = "";
    $("inp-repeat").disabled = false;
  } else {
    $("inp-repeat").value = 1;
    row.style.display = "none";
  }
}

/* ---------- History: snapshot helpers ------------------------------------- */
function takeSnapshot() {
  return {
    jobId:   state.jobId,
    built:   state.built,
    xyzStr:  state.xyzStr,
    box:     state.box ? state.box.slice() : null,
    kind:    state.kind,
    summary: state.summary ? { ...state.summary } : null,
    walls:   state.walls,
    warning: state.warning,
    torsionApplied: state.torsionApplied,
    atomsBase: state.atomsBase,
    deformDesc: state.deformDesc,
    repsValue: parseInt($("inp-repeat").value, 10) || 1,
    repsHidden: $("reps-row").style.display === "none",
  };
}
function applySnapshot(s) {
  state.jobId   = s.jobId;
  state.built   = s.built;
  state.xyzStr  = s.xyzStr;
  state.box     = s.box;
  state.kind    = s.kind;
  state.summary = s.summary;
  state.walls   = s.walls;
  state.warning = s.warning;
  state.torsionApplied = s.torsionApplied;
  state.atomsBase  = s.atomsBase;
  state.deformDesc = s.deformDesc;
  $("inp-repeat").value = s.repsValue ?? 1;
  setRepsVisible(!s.repsHidden);
  // Restore viewer ×z lock status to match the snapshot.
  const vrz = $("view-rep-z");
  if (s.torsionApplied) {
    vrz.value = 1; vrz.disabled = true;
    vrz.title = t("viewer.repz.locked");
  } else {
    vrz.disabled = false;
    vrz.title = t("viewer.repz.title");
  }
  updateAfterBuild(/*pushHistory=*/false);
  if (state.xyzStr) loadViewer(state.xyzStr, state.box);
  else { clearViewer(); }
}
function pushUndo() {
  state.undoStack.push(takeSnapshot());
  if (state.undoStack.length > 20) state.undoStack.shift();
  state.redoStack.length = 0;
  refreshHistoryButtons();
}
function refreshHistoryButtons() {
  $("undo-btn").disabled = state.undoStack.length === 0;
  $("redo-btn").disabled = state.redoStack.length === 0;
}
$("undo-btn").addEventListener("click", () => {
  if (!state.undoStack.length) return;
  const prev = state.undoStack.pop();
  state.redoStack.push(takeSnapshot());
  applySnapshot(prev);
  refreshHistoryButtons();
  toast(t("t.undone"), "info", 1800);
});
$("redo-btn").addEventListener("click", () => {
  if (!state.redoStack.length) return;
  const next = state.redoStack.pop();
  state.undoStack.push(takeSnapshot());
  applySnapshot(next);
  refreshHistoryButtons();
  toast(t("t.redone"), "info", 1800);
});
document.addEventListener("keydown", e => {
  if (e.target.matches("input, textarea, select")) return;
  const ctrl = e.ctrlKey || e.metaKey;
  if (ctrl && !e.shiftKey && e.key.toLowerCase() === "z") {
    e.preventDefault(); $("undo-btn").click();
  } else if (ctrl && (e.key.toLowerCase() === "y" || (e.shiftKey && e.key.toLowerCase() === "z"))) {
    e.preventDefault(); $("redo-btn").click();
  }
});

/* ---------- Chain-warning prompts ----------------------------------------- */
const CHAIN_WARNING_KEYS = {
  mwnt:   { bundle: "chain.mwnt.bundle",  mwnt: "chain.mwnt.mwnt",  torsion: "chain.mwnt.torsion" },
  bundle: { torsion: "chain.bundle.torsion" },
  axial:  { torsion: "chain.axial.torsion" },
  torsion:{ bundle:  "chain.torsion.bundle" },
};
function confirmChain(nextOp) {
  const cur = state.kind;
  const key = (CHAIN_WARNING_KEYS[nextOp] || {})[cur];
  if (!key) return Promise.resolve(true);
  return new Promise(resolve => {
    $("confirm-title").textContent = t("chain.title", { op: nextOp });
    $("confirm-text").innerHTML    = t(key);
    $("confirm-overlay").classList.add("open");
    const yes = $("confirm-yes"), no = $("confirm-no");
    function cleanup() {
      $("confirm-overlay").classList.remove("open");
      yes.removeEventListener("click", okFn);
      no.removeEventListener("click", noFn);
    }
    function okFn() { cleanup(); resolve(true); }
    function noFn() { cleanup(); resolve(false); }
    yes.addEventListener("click", okFn);
    no.addEventListener("click", noFn);
  });
}

/* ---------- Build (SWNT) --------------------------------------------------- */
$("build-btn").addEventListener("click", runBuild);

async function runBuild() {
  const n = parseInt($("inp-n").value, 10), m = parseInt($("inp-m").value, 10);
  if (n === 0 && m === 0) { toast(t("t.degen"), "err"); return; }

  pushUndo();          // permitir desfazer um build novo
  const btn = $("build-btn");
  btn.disabled = true;
  btn.innerHTML = `<span class="spinner" style="width:14px;height:14px;margin-right:8px;"></span>${t("step3.building")}`;

  const body = {
    n, m,
    n_repeat:    parseInt($("inp-repeat").value, 10) || 1,
    vacuum:      parseFloat($("inp-vacuum").value) || 10,
    roll_inward: $("inp-rollin").checked,
    ...(Object.keys(state.bondCutoffs).length ? { bond_cutoffs: state.bondCutoffs } : {}),
    ...(state.fileId  ? { file_id: state.fileId  } : {}),
    ...(state.example ? { example: state.example } : {}),
  };
  try {
    const res = await apiJSON("POST", "api/build", body);
    state.jobId   = res.job_id;
    state.xyzStr  = res.xyz;
    state.box     = res.box || null;
    state.built   = true;
    state.kind    = "swnt";
    state.summary = res;
    state.walls   = null;
    state.warning = res.warning || null;
    state.torsionApplied = false;
    state.deformDesc = null;
    state.repsAtBuild = res.n_repeat || 1;
    state.atomsBase = res.n_atoms / (res.n_repeat || 1);
    setRepsVisible(true);
    // Re-enable the viewer's ×z replication: a fresh build restores
    // axial periodicity even if a previous torsion had locked it.
    const vrz = $("view-rep-z");
    vrz.disabled = false;
    vrz.title = t("viewer.repz.title");
    updateAfterBuild();
    if (res.xyz) loadViewer(res.xyz, res.box);
    toast(`${t("t.buildOK.swnt")}: ${res.n_atoms.toLocaleString()} ${t("step4.atoms").toLowerCase()}`, "ok");
  } catch (err) {
    toast(t("t.buildFail") + ": " + err.message, "err", 5000);
    // Reverte o snapshot que acabamos de empilhar — nada mudou de fato
    state.undoStack.pop();
    refreshHistoryButtons();
  } finally {
    btn.disabled = false;
    btn.textContent = t("step3.build");
  }
}

function updateAfterBuild(pushHistory = true) {
  const s = state.summary;
  if (s && s.n !== undefined) {
    $("sg-nm").textContent     = `(${s.n}, ${s.m})`;
    $("sg-d").textContent      = s.diameter !== undefined ? `${s.diameter} Å` : "—";
    $("sg-l").textContent      = s.length   !== undefined ? `${(+s.length).toFixed(2)} Å` : "—";
    $("sg-nat").textContent    = (s.n_atoms || 0).toLocaleString();
    $("sg-strain").textContent = s.strain   !== undefined ? `${s.strain.toFixed(4)} %` : "—";
    $("sg-theta").textContent  = s.theta_deg !== undefined ? `${s.theta_deg}°` : "—";
  } else {
    ["sg-nm", "sg-d", "sg-l", "sg-nat", "sg-strain", "sg-theta"]
      .forEach(id => $(id).textContent = "—");
  }
  const w = $("warning-box");
  if (state.warning) {
    w.textContent = "⚠ " + state.warning;
    w.classList.add("visible");
  } else {
    w.classList.remove("visible");
    w.textContent = "";
  }
  // Deformação aplicada — mostra abaixo do summary se houver
  if (state.deformDesc && state.deformDesc !== "none") {
    w.textContent = `${t("t.deformChip")}: ${state.deformDesc}`
                 + (state.warning ? "\n⚠ " + state.warning : "");
    w.classList.add("visible");
  }
  renderWallPlan();
  updateActionButtons();
  updateRepsInfo();
  refreshHistoryButtons();
}

function renderWallPlan() {
  const el = $("wall-plan");
  if (!state.walls || !state.walls.length) {
    el.classList.add("hidden"); el.textContent = "";
    return;
  }
  const lines = [];
  lines.push(`${state.walls.length} ${t("t.wallsLabel")}`);
  lines.push("");
  lines.push(" idx   k   (n, m)        D (Å)    gap (Å)");
  state.walls.forEach(w => {
    const gap = (w.actual_spacing == null) ? "  —  " : w.actual_spacing.toFixed(3);
    lines.push(
      ` ${String(w.index).padStart(3)}` +
      ` ${String(w.k).padStart(4)}` +
      `  (${w.n},${w.m})`.padEnd(14) +
      ` ${w.diameter.toFixed(3).padStart(8)}` +
      ` ${gap.padStart(9)}`,
    );
  });
  el.innerHTML = `<div class="wp-title">${t("t.wallPlanTitle")}</div>${lines.join("\n")}`;
  el.classList.remove("hidden");
}

function updateActionButtons() {
  const built = state.built;
  qsa("#download-buttons .dl-btn").forEach(b => b.disabled = !built);
  $("download-open-btn").disabled = !built;
  $("mwnt-open-btn").disabled     = !built;
  $("bundle-open-btn").disabled   = !built;
  $("deform-open-btn").disabled   = !built;
  $("analysis-open-btn").disabled = !built;
  $("methods-open-btn").disabled  = !built;
  $("dft-open-btn").disabled      = !built;
}

/* ---------- 3D viewer ------------------------------------------------------ */
function replicateXYZ(xyzStr, nRep, boxZ) {
  if (nRep <= 1 || !boxZ) return xyzStr;
  const lines = xyzStr.trim().split("\n");
  const nAtoms = parseInt(lines[0], 10);
  const comment = lines[1];
  const atomLines = lines.slice(2);
  const newAtoms = [];
  for (let i = 0; i < nRep; i++) {
    atomLines.forEach(line => {
      const p = line.trim().split(/\s+/);
      newAtoms.push(`${p[0]} ${p[1]} ${p[2]} ${(parseFloat(p[3]) + i * boxZ).toFixed(6)}`);
    });
  }
  return `${nAtoms * nRep}\n${comment}\n${newAtoms.join("\n")}`;
}

function loadViewer(xyzStr, box) {
  $("viewer-empty").style.display = "none";
  const container = $("viewer-3d");
  if (state.viewer) { try { state.viewer.clear(); } catch (_) {} }

  // Defensive: after torsion, axial periodicity is broken and the backend
  // already includes the chosen reps in the returned XYZ.  Replicating again
  // in JS would (a) hide the twist visually (identical translated copies)
  // and (b) blow up the atom count quadratically when the user had left
  // ``view-rep-z`` at a large value.  Force nRep=1 in that case regardless
  // of the spinbox value.
  const vrzVal = parseInt($("view-rep-z").value, 10) || 1;
  const nRep   = state.torsionApplied ? 1 : vrzVal;
  const boxZ   = box ? box[2] : null;
  const xyz    = replicateXYZ(xyzStr, nRep, boxZ);

  const viewer = $3Dmol.createViewer(container, {
    backgroundColor: getViewerBg(), antialias: true,
  });
  state.viewer = viewer;
  viewer.addModel(xyz, "xyz");
  applyViewStyle(viewer, state.viewStyle);
  if (box && state.showBox) {
    const isDark = document.documentElement.getAttribute("data-theme") === "dark";
    drawSimBox(viewer, box[0], box[1], box[2] * nRep, isDark);
  }
  viewer.zoomTo(); viewer.render();
  $("viewer-controls").classList.remove("hidden");
  $("viewer-controls").style.display = "flex";
}

function clearViewer() {
  if (state.viewer) { try { state.viewer.clear(); state.viewer.render(); } catch (_) {} }
  $("viewer-empty").style.display = "flex";
  $("viewer-controls").classList.add("hidden");
  $("viewer-controls").style.display = "none";
}

function drawSimBox(viewer, Lx, Ly, Lz, dark) {
  const color = dark ? "#CBD5E1" : "#475569";
  const C = [
    [0,0,0],[Lx,0,0],[Lx,Ly,0],[0,Ly,0],
    [0,0,Lz],[Lx,0,Lz],[Lx,Ly,Lz],[0,Ly,Lz],
  ];
  const edges = [[0,1],[1,2],[2,3],[3,0],[4,5],[5,6],[6,7],[7,4],[0,4],[1,5],[2,6],[3,7]];
  edges.forEach(([a, b]) => {
    viewer.addCylinder({
      start: { x: C[a][0], y: C[a][1], z: C[a][2] },
      end:   { x: C[b][0], y: C[b][1], z: C[b][2] },
      radius: 0.06, color, opacity: 0.6, fromCap: false, toCap: false,
    });
  });
}
function applyViewStyle(viewer, styleIdx) {
  viewer.setStyle({}, {});
  if      (styleIdx === 0) viewer.setStyle({}, { sphere: { radius: 0.25 }, stick: { radius: 0.12 } });
  else if (styleIdx === 1) viewer.setStyle({}, { sphere: {} });
  else                     viewer.setStyle({}, { stick: { radius: 0.15 } });
  viewer.render();
}

$("view-reset-btn").addEventListener("click", () => {
  if (state.viewer) { state.viewer.zoomTo(); state.viewer.render(); }
});
$("view-style-btn").addEventListener("click", () => {
  if (!state.viewer) return;
  state.viewStyle = (state.viewStyle + 1) % 3;
  applyViewStyle(state.viewer, state.viewStyle);
});
$("view-bg-btn").addEventListener("click", () => {
  if (!state.viewer) return;
  const dark = document.documentElement.getAttribute("data-theme") === "dark";
  state.viewBg = dark ? "#F8FAFC" : "#060B17";
  state.viewer.setBackgroundColor(state.viewBg);
  state.viewer.render();
});
$("view-box-btn").addEventListener("click", () => {
  state.showBox = !state.showBox;
  $("view-box-btn").style.borderColor = state.showBox ? "var(--c-brand)" : "";
  if (state.xyzStr) loadViewer(state.xyzStr, state.box);
});
["input", "change"].forEach(ev =>
  $("view-rep-z").addEventListener(ev, () => {
    if (state.xyzStr) loadViewer(state.xyzStr, state.box);
  })
);

/* ---------- Downloads ----------------------------------------------------- */
// Os 9 botões originais ficam no DOM (escondidos via .hidden no container) só
// para preservar a lógica de habilitar/desabilitar centralizada em
// updateActionButtons().  O usuário acessa-os pela janela modal "Baixar
// estrutura", que delega cada clique para o botão correspondente.
function downloadFmt(fmt) {
  if (!state.jobId) { toast(t("t.buildFirst"), "err"); return; }
  const a = document.createElement("a");
  a.href = `api/download/${state.jobId}/${fmt}`;
  a.download = ""; a.click();
}
qsa("#download-buttons .dl-btn").forEach(b => {
  b.addEventListener("click", () => downloadFmt(b.dataset.fmt));
});
qsa('[data-modal-fmt]').forEach(b => {
  b.addEventListener("click", () => downloadFmt(b.dataset.modalFmt));
});

$("download-open-btn").addEventListener("click", () => {
  if (!state.built) { toast(t("t.buildFirst"), "err"); return; }
  $("download-overlay").classList.add("open");
});
$("download-close-btn").addEventListener("click",
  () => $("download-overlay").classList.remove("open"));

/* ---------- Generic advanced op runner ------------------------------------ */
async function runAdvancedOp(endpoint, body, sourceLabelKey, opts = {}) {
  const { chainable = true, beforeRequest = null, onSuccess = null } = opts;
  if (!state.fileId && !state.example) {
    toast(t("t.loadStructFirst"), "err");
    return null;
  }
  if (!state.built) {
    toast(t("t.buildFirst"), "err");
    return null;
  }
  const sourceLabel = sourceLabelKey;
  if (beforeRequest) {
    const proceed = await beforeRequest();
    if (!proceed) return null;
  }
  pushUndo();   // snapshot before mutation
  const payload = {
    ...body,
    n: parseInt($("inp-n").value, 10),
    m: parseInt($("inp-m").value, 10),
    ...(state.fileId  ? { file_id: state.fileId  } : {}),
    ...(state.example ? { example: state.example } : {}),
    ...(chainable && state.jobId ? { from_job_id: state.jobId } : {}),
  };
  try {
    const res = await apiJSON("POST", endpoint, payload);
    if (res.job_id) {
      state.jobId = res.job_id;
      state.box   = res.box || null;
      try {
        const xyzResp = await api("GET", `api/xyz/${res.job_id}`);
        state.xyzStr = await xyzResp.text();
      } catch (_) { state.xyzStr = null; }
      state.built = true;
    }
    // IMPORTANT: onSuccess must run BEFORE loadViewer.
    //
    // The deform callback flips `state.torsionApplied = true` and forces the
    // viewer's ×z replication spinbox to 1 — that prevents loadViewer (which
    // re-replicates the XYZ in pure JS via ``replicateXYZ``) from compounding
    // the axial replication that the backend has already baked into the
    // returned coordinates.  Running loadViewer before onSuccess caused the
    // returned 30×-replicated supercell to be replicated *another* 30× on the
    // client, blowing the atom count up by 900×: a hex7 bundle (~1 k atoms)
    // turned into ~10⁶ atoms processed by a synchronous string-split loop,
    // hanging the browser tab for minutes.
    if (onSuccess) onSuccess(res);
    // Mirror the resulting torsion state onto the ×z control: any op that
    // isn't a torsion gives back axial periodicity, so the spinbox unlocks.
    if (!state.torsionApplied) {
      const vrz = $("view-rep-z");
      vrz.disabled = false;
      vrz.title = t("viewer.repz.title");
    }
    if (state.xyzStr) loadViewer(state.xyzStr, state.box);
    toast(t("t.opOK", { label: sourceLabel, n: (res.n_atoms || 0).toLocaleString() }), "ok");
    updateAfterBuild(/*pushHistory=*/false);
    return res;
  } catch (err) {
    toast(t("t.opFail", { label: sourceLabel }) + ": " + err.message, "err", 6000);
    state.undoStack.pop();
    refreshHistoryButtons();
    return null;
  }
}

/* ---------- MWNT ----------------------------------------------------------- */
$("mwnt-open-btn").addEventListener("click", () => {
  $("mwnt-plan-modal").classList.add("hidden");
  $("mwnt-plan-modal").textContent = "";
  $("mwnt-overlay").classList.add("open");
});
$("mwnt-cancel-btn").addEventListener("click",
  () => $("mwnt-overlay").classList.remove("open"));

$("mwnt-build-btn").addEventListener("click", async () => {
  const btn = $("mwnt-build-btn");
  btn.disabled = true;
  btn.innerHTML = `<span class="spinner" style="width:14px;height:14px;margin-right:8px;"></span>${t("mwnt.going")}`;
  try {
    const res = await runAdvancedOp("api/mwnt", {
      n_walls:            parseInt($("mwnt-walls").value, 10),
      interlayer_spacing: parseFloat($("mwnt-spacing").value),
      vacuum:             parseFloat($("mwnt-vacuum").value),
      roll_inward:        $("mwnt-rollin").checked,
    }, "MWNT", {
      chainable: false,
      beforeRequest: () => confirmChain("mwnt"),
      onSuccess: (r) => {
        state.kind  = "mwnt";
        state.walls = r.walls;
        state.warning = null;
        state.torsionApplied = false;
        state.atomsBase = r.n_atoms;
        state.summary = {
          n: parseInt($("inp-n").value, 10),
          m: parseInt($("inp-m").value, 10),
          n_atoms: r.n_atoms,
          length: r.box ? r.box[2] : 0,
        };
        // Mostra o plano de paredes também dentro do modal
        if (r.walls && r.walls.length) {
          const lines = [
            `${r.walls.length} ${t("t.wallsLabel")} — ${r.requested_spacing.toFixed(2)} Å`
              + `, ⌀ ${r.mean_spacing.toFixed(2)} Å`,
            "",
            " idx   k   (n, m)        D (Å)    gap (Å)",
          ];
          r.walls.forEach(w => {
            const gap = (w.actual_spacing == null) ? "  —  " : w.actual_spacing.toFixed(3);
            lines.push(
              ` ${String(w.index).padStart(3)}` +
              ` ${String(w.k).padStart(4)}` +
              `  (${w.n},${w.m})`.padEnd(14) +
              ` ${w.diameter.toFixed(3).padStart(8)}` +
              ` ${gap.padStart(9)}`);
          });
          $("mwnt-plan-modal").textContent = lines.join("\n");
          $("mwnt-plan-modal").classList.remove("hidden");
        }
      },
    });
    if (res) $("mwnt-overlay").classList.remove("open");
  } finally {
    btn.disabled = false;
    btn.textContent = t("mwnt.go");
  }
});

/* ---------- Bundle --------------------------------------------------------- */
$("bundle-open-btn").addEventListener("click",
  () => $("bundle-overlay").classList.add("open"));
$("bundle-cancel-btn").addEventListener("click",
  () => $("bundle-overlay").classList.remove("open"));
$("bundle-geometry").addEventListener("change", () => {
  $("bundle-grid-row").style.display =
    ($("bundle-geometry").value === "grid") ? "" : "none";
});
$("bundle-build-btn").addEventListener("click", async () => {
  const btn = $("bundle-build-btn");
  btn.disabled = true;
  btn.innerHTML = `<span class="spinner" style="width:14px;height:14px;margin-right:8px;"></span>${t("bundle.going")}`;
  try {
    const res = await runAdvancedOp("api/bundle", {
      geometry: $("bundle-geometry").value,
      spacing:  parseFloat($("bundle-spacing").value),
      vacuum:   parseFloat($("bundle-vacuum").value),
      nx:       parseInt($("bundle-nx").value, 10),
      ny:       parseInt($("bundle-ny").value, 10),
      n_repeat: parseInt($("bundle-repeat").value, 10),
    }, "Bundle", {
      beforeRequest: () => confirmChain("bundle"),
      onSuccess: (r) => {
        state.kind = "bundle";
        state.walls = null;
        state.warning = null;
        state.atomsBase = r.n_atoms;
        state.summary = {
          n: parseInt($("inp-n").value, 10),
          m: parseInt($("inp-m").value, 10),
          n_atoms: r.n_atoms,
          length:  r.box ? r.box[2] : 0,
        };
      },
    });
    if (res) $("bundle-overlay").classList.remove("open");
  } finally {
    btn.disabled = false;
    btn.textContent = t("bundle.go");
  }
});

/* ---------- Deform / Torsion ----------------------------------------------- */
// Limite de átomos do backend para /api/deform — mantém em sync com main.py
const DEFORM_MAX_ATOMS = 250_000;

function updateDeformProjection() {
  const reps = Math.max(1, parseInt($("deform-repeat").value, 10) || 1);
  const base = state.atomsBase || (state.summary && state.summary.n_atoms) || 0;
  const proj = base * reps;
  const el = $("deform-projection");
  if (!base) {
    el.classList.add("hidden");
    return;
  }
  el.classList.remove("hidden");
  el.classList.toggle("warning", proj > DEFORM_MAX_ATOMS);
  const tag = proj > DEFORM_MAX_ATOMS ? "⚠" : "ℹ";
  el.textContent =
    `${tag} ${t("t.projInfo")}: ${reps} × ${base.toLocaleString()} = ${proj.toLocaleString()}`
    + (proj > DEFORM_MAX_ATOMS
        ? `\n${t("t.projOver", { max: DEFORM_MAX_ATOMS.toLocaleString() })}`
        : "");
}

$("deform-open-btn").addEventListener("click", () => {
  $("deform-warning").classList.add("hidden");
  $("deform-warning").classList.remove("warning");
  $("deform-warning").textContent = "";

  // The torsion is always applied to the supercell that the user is currently
  // visualising (the ×z spinbox in the viewer header).  We mirror that value
  // into the hidden ``deform-repeat`` field so the rest of the pipeline
  // (Apply handler, undo/redo snapshot, projection check) stays unchanged.
  // The field used to be editable in the modal — that caused user confusion
  // (looked like an extra replication on top of what was already on screen)
  // so the spinbox was removed; the value is read from the viewer instead.
  const viewerRep = parseInt($("view-rep-z").value, 10) || 1;
  $("deform-repeat").value = Math.min(999, viewerRep);
  updateDeformProjection();
  $("deform-overlay").classList.add("open");
});
$("deform-cancel-btn").addEventListener("click",
  () => $("deform-overlay").classList.remove("open"));

$("deform-build-btn").addEventListener("click", async () => {
  const strainPct = parseFloat($("deform-strain").value) || 0;
  const radialPct = parseFloat($("deform-radial").value) || 0;
  const twist     = parseFloat($("deform-twist").value)  || 0;
  const reps      = parseInt($("deform-repeat").value, 10) || 1;

  // Bloqueio preventivo: não envia requisição se a projeção já excede o cap
  const base = state.atomsBase || (state.summary && state.summary.n_atoms) || 0;
  if (base && reps * base > DEFORM_MAX_ATOMS) {
    toast(t("t.deformProjExceeds", {
      n: (reps * base).toLocaleString(),
      max: DEFORM_MAX_ATOMS.toLocaleString(),
    }), "err", 7000);
    return;
  }

  // Chain warning depende de qual deformação foi escolhida
  const nextOp = (Math.abs(twist) > 1e-9) ? "torsion"
              : (Math.abs(strainPct) > 1e-9 ? "axial" : "radial");

  const btn = $("deform-build-btn");
  btn.disabled = true;
  btn.innerHTML = `<span class="spinner" style="width:14px;height:14px;margin-right:8px;"></span>${t("deform.going")}`;
  try {
    const res = await runAdvancedOp("api/deform", {
      axial_strain:  strainPct / 100.0,
      twist_rate:    twist,
      radial_strain: radialPct / 100.0,
      vacuum:        parseFloat($("deform-vacuum").value),
      z_vacuum:      parseFloat($("deform-zvac").value),
      n_repeat:      parseInt($("deform-repeat").value, 10),
    }, "Deformation", {
      beforeRequest: () => confirmChain(nextOp),
      onSuccess: (r) => {
        state.kind = (Math.abs(twist) > 1e-9) ? "torsion"
                   : (Math.abs(strainPct) > 1e-9 ? "axial" : "radial");
        state.torsionApplied = (Math.abs(twist) > 1e-9);
        state.walls   = null;
        state.warning = r.warning || null;
        state.deformDesc = r.description || null;
        state.atomsBase = r.n_atoms;
        state.summary = {
          n: parseInt($("inp-n").value, 10),
          m: parseInt($("inp-m").value, 10),
          n_atoms: r.n_atoms,
          length:  r.box ? r.box[2] : 0,
        };
        // Torsion breaks axial periodicity. Hide the sidebar Reps row and
        // also clamp the viewer's ×z replication to 1 — replicating a
        // twisted slice translationally yields N identical copies (no
        // continuation of the twist) and gives the false impression that
        // torsion was never applied.
        if (state.torsionApplied) {
          setRepsVisible(false);
          const vrz = $("view-rep-z");
          vrz.value = 1; vrz.disabled = true;
          vrz.title = t("viewer.repz.locked");
        }
        const msgs = [];
        if (r.description && r.description !== "none") msgs.push(`${t("t.deformApplied")}: ${r.description}`);
        if (r.warning) msgs.push(r.warning);
        if (msgs.length) {
          $("deform-warning").textContent = msgs.join("\n\n");
          $("deform-warning").classList.remove("hidden");
          if (r.warning) $("deform-warning").classList.add("warning");
        }
      },
    });
    if (res) $("deform-overlay").classList.remove("open");
  } finally {
    btn.disabled = false;
    btn.textContent = t("deform.go");
  }
});

/* ---------- Analysis (bond histogram + electronic + symmetry) -------------- */
$("analysis-open-btn").addEventListener("click",
  () => $("analysis-overlay").classList.add("open"));
$("analysis-close-btn").addEventListener("click",
  () => $("analysis-overlay").classList.remove("open"));

$("analysis-run-btn").addEventListener("click", async () => {
  if (!state.built) { toast(t("t.buildFirst"), "err"); return; }
  const n = parseInt($("inp-n").value, 10), m = parseInt($("inp-m").value, 10);
  const body = {
    n, m,
    vacuum:      parseFloat($("inp-vacuum").value) || 10,
    bond_cutoff: parseFloat($("analysis-cutoff").value) || 2.0,
    ...(state.fileId  ? { file_id: state.fileId  } : {}),
    ...(state.example ? { example: state.example } : {}),
  };
  const btn = $("analysis-run-btn");
  btn.disabled = true;
  btn.innerHTML = `<span class="spinner" style="width:13px;height:13px;margin-right:6px;"></span>${t("analysis.running")}`;
  try {
    const res = await apiJSON("POST", "api/analysis", body);
    renderAnalysis(res);
  } catch (err) {
    toast(t("t.analysisFail") + ": " + err.message, "err", 6000);
  } finally {
    btn.disabled = false;
    btn.textContent = t("analysis.run");
  }
});

function renderAnalysis(res) {
  const ba = res.bond_analysis;
  const elec = res.electronic_label;
  const sym  = res.symmetry;
  const html = `
    <b>${t("analysis.bondStats")}</b> — N=${ba.n_bonds.toLocaleString()} ${t("analysis.bonds")}
    · ${t("analysis.mean")}=${ba.mean.toFixed(4)} Å · ${t("analysis.std")}=${ba.std.toFixed(4)} Å
    · ${t("analysis.min")}=${ba.min.toFixed(3)} Å · ${t("analysis.max")}=${ba.max.toFixed(3)} Å<br>
    <b>${t("analysis.pairs")}</b> ${ba.species.join(", ") || "—"}<br>
    <b>${t("analysis.elec")}</b> ${elec}<br>
    <b>${t("analysis.sym")}</b> ${sym.description}
  `;
  const sumEl = $("analysis-summary");
  sumEl.innerHTML = html;
  sumEl.classList.remove("hidden");

  // Plotly histogram (per-pair stacked if multiple species pairs)
  const c = plotColors();
  const pairs = Array.from(new Set(ba.pairs));
  const traces = pairs.map(p => ({
    type: "histogram",
    x: ba.distances.filter((_, i) => ba.pairs[i] === p),
    name: p, opacity: 0.75,
    xbins: { size: 0.02 },
  }));
  const layout = {
    paper_bgcolor: "transparent", plot_bgcolor: "transparent",
    margin: { l: 50, r: 20, t: 10, b: 40 },
    barmode: pairs.length > 1 ? "stack" : "overlay",
    xaxis: { title: { text: t("analysis.xLabel"), font: { color: c.text } },
             gridcolor: c.grid, color: c.text, tickfont: { color: c.text } },
    yaxis: { title: { text: t("analysis.yLabel"), font: { color: c.text } },
             gridcolor: c.grid, color: c.text, tickfont: { color: c.text } },
    legend: { orientation: "h", x: 0, y: 1.12, font: { size: 10, color: c.text } },
    font: { family: "Inter, sans-serif", color: c.text },
  };
  Plotly.react("analysis-plot", traces, layout, { responsive: true, displayModeBar: false });
}

/* ---------- Methods text --------------------------------------------------- */
$("methods-open-btn").addEventListener("click",
  () => $("methods-overlay").classList.add("open"));
$("methods-close-btn").addEventListener("click",
  () => $("methods-overlay").classList.remove("open"));

$("methods-run-btn").addEventListener("click", async () => {
  if (!state.built) { toast(t("t.buildFirst"), "err"); return; }
  const body = {
    n: parseInt($("inp-n").value, 10),
    m: parseInt($("inp-m").value, 10),
    vacuum:   parseFloat($("inp-vacuum").value) || 10,
    cite_key: $("methods-citekey").value || "Pereira2026",
    deform_desc: state.deformDesc || "",
    n_walls: (state.walls && state.walls.length) ? state.walls.length : 1,
    wall_info: state.walls
      ? state.walls.map(w => `wall ${w.index}: (${w.n},${w.m}) D=${w.diameter.toFixed(2)} Å`).join("; ")
      : "",
    ...(state.fileId   ? { file_id:     state.fileId  } : {}),
    ...(state.example  ? { example:     state.example } : {}),
    ...(state.jobId    ? { from_job_id: state.jobId   } : {}),
  };
  const btn = $("methods-run-btn");
  btn.disabled = true;
  btn.innerHTML = `<span class="spinner" style="width:13px;height:13px;margin-right:6px;"></span>${t("methods.going")}`;
  try {
    const res = await apiJSON("POST", "api/methods", body);
    const out = $("methods-output");
    out.textContent = res.text || "";
    out.classList.remove("hidden");
    $("methods-copy-btn").disabled = false;
  } catch (err) {
    toast(t("t.methodsFail") + ": " + err.message, "err", 6000);
  } finally {
    btn.disabled = false;
    btn.textContent = t("methods.gen");
  }
});

$("methods-copy-btn").addEventListener("click", () => {
  const text = $("methods-output").textContent;
  if (!text) return;
  navigator.clipboard.writeText(text)
    .then(() => toast(t("methods.copied"), "ok"))
    .catch(() => toast(t("t.copyFail"), "err"));
});

/* ---------- DFT inputs ----------------------------------------------------- */
$("dft-open-btn").addEventListener("click",
  () => $("dft-overlay").classList.add("open"));
$("dft-close-btn").addEventListener("click",
  () => $("dft-overlay").classList.remove("open"));

let _dftFiles = null;
let _dftActive = null;

$("dft-run-btn").addEventListener("click", async () => {
  if (!state.built) { toast(t("t.buildFirst"), "err"); return; }
  const body = {
    n: parseInt($("inp-n").value, 10),
    m: parseInt($("inp-m").value, 10),
    vacuum: parseFloat($("inp-vacuum").value) || 10,
    code:   $("dft-code").value,
    ...(state.fileId   ? { file_id:     state.fileId  } : {}),
    ...(state.example  ? { example:     state.example } : {}),
    ...(state.jobId    ? { from_job_id: state.jobId   } : {}),
  };
  const btn = $("dft-run-btn");
  btn.disabled = true;
  btn.innerHTML = `<span class="spinner" style="width:13px;height:13px;margin-right:6px;"></span>${t("dft.going")}`;
  try {
    const res = await apiJSON("POST", "api/dft_inputs", body);
    _dftFiles = res.files || {};
    renderDFTTabs(_dftFiles);
    $("dft-copy-btn").disabled = false;
  } catch (err) {
    toast(t("t.dftFail") + ": " + err.message, "err", 6000);
  } finally {
    btn.disabled = false;
    btn.textContent = t("dft.gen");
  }
});

function renderDFTTabs(files) {
  const tabs = $("dft-tabs");
  const out  = $("dft-output");
  tabs.innerHTML = "";
  const names = Object.keys(files);
  if (!names.length) {
    tabs.classList.add("hidden");
    out.classList.add("hidden");
    return;
  }
  _dftActive = names[0];
  names.forEach(name => {
    const b = document.createElement("button");
    b.textContent = name;
    if (name === _dftActive) b.classList.add("active");
    b.addEventListener("click", () => {
      _dftActive = name;
      tabs.querySelectorAll("button").forEach(x => x.classList.remove("active"));
      b.classList.add("active");
      out.textContent = files[name];
    });
    tabs.appendChild(b);
  });
  tabs.classList.remove("hidden");
  out.textContent = files[_dftActive];
  out.classList.remove("hidden");
}

$("dft-copy-btn").addEventListener("click", () => {
  if (!_dftFiles || !_dftActive) return;
  navigator.clipboard.writeText(_dftFiles[_dftActive])
    .then(() => toast(t("dft.copied", { name: _dftActive }), "ok"))
    .catch(() => toast(t("t.copyFail"), "err"));
});

/* ---------- Bond cutoffs --------------------------------------------------- */
$("btn-cutoffs").addEventListener("click", () => {
  renderCutoffs();
  $("cutoffs-overlay").classList.add("open");
});
$("cutoffs-close-btn").addEventListener("click",
  () => $("cutoffs-overlay").classList.remove("open"));

function renderCutoffs() {
  const list = $("cutoffs-list");
  list.innerHTML = "";
  Object.entries(state.bondCutoffs).forEach(([pair, val]) => {
    const row = document.createElement("div");
    row.className = "cutoff-row";
    row.innerHTML =
      `<span>${pair}</span>`
      + `<input type="number" step="0.05" min="0.5" max="5" value="${val}" />`
      + `<button title="Remove">✕</button>`;
    row.querySelector("input").addEventListener("change", e => {
      state.bondCutoffs[pair] = parseFloat(e.target.value);
    });
    row.querySelector("button").addEventListener("click", () => {
      delete state.bondCutoffs[pair];
      renderCutoffs();
    });
    list.appendChild(row);
  });
  if (!Object.keys(state.bondCutoffs).length) {
    list.innerHTML = `<div style="font-size:.78rem;color:var(--c-muted);">${t("cutoffs.empty")}</div>`;
  }
}

$("cutoff-add-btn").addEventListener("click", () => {
  const pair = ($("cutoff-new-pair").value || "").trim();
  const val  = parseFloat($("cutoff-new-val").value);
  if (!pair || !Number.isFinite(val)) {
    toast(t("cutoffs.invalid"), "err");
    return;
  }
  state.bondCutoffs[pair] = val;
  $("cutoff-new-pair").value = "";
  $("cutoff-new-val").value  = "";
  renderCutoffs();
});

/* ---------- Batch ---------------------------------------------------------- */
$("batch-open-btn").addEventListener("click",
  () => $("batch-overlay").classList.add("open"));
$("batch-cancel-btn").addEventListener("click",
  () => $("batch-overlay").classList.remove("open"));
qsa('input[name="batch-type"]').forEach(r => {
  r.addEventListener("change", () => {
    $("batch-range-row").style.display = r.value === "all" ? "none" : "";
  });
});
$("batch-build-btn").addEventListener("click", runBatch);

async function runBatch() {
  if (!state.fileId && !state.example) { toast(t("t.loadStructFirst"), "err"); return; }
  const batchType = document.querySelector('input[name="batch-type"]:checked').value;
  let chiralities = [];
  if (batchType === "all") {
    if (!state.polarData) { toast(t("t.batchNeedMap"), "err"); return; }
    chiralities = getFilteredPoints(state.polarData).map(p => [p.n, p.m]);
  } else {
    const nFrom = parseInt($("batch-n-from").value, 10);
    const nTo   = parseInt($("batch-n-to").value, 10);
    if (nFrom > nTo) { toast(t("t.batchOrder"), "err"); return; }
    for (let n = nFrom; n <= nTo; n++)
      chiralities.push(batchType === "armchair" ? [n, n] : [n, 0]);
  }
  if (!chiralities.length) { toast(t("t.batchEmpty"), "err"); return; }
  if (chiralities.length > 100) { toast(t("t.batchMax"), "err"); return; }

  const btn = $("batch-build-btn");
  btn.disabled = true;
  btn.innerHTML = `<span class="spinner" style="width:14px;height:14px;margin-right:8px;"></span>${t("batch.going", { n: chiralities.length })}`;
  const body = {
    chiralities,
    n_repeat:    parseInt($("batch-repeat").value, 10),
    vacuum:      parseFloat($("batch-vacuum").value),
    roll_inward: $("batch-rollin").checked,
    ...(state.fileId  ? { file_id: state.fileId  } : {}),
    ...(state.example ? { example: state.example } : {}),
  };
  try {
    const res  = await api("POST", "api/batch", body);
    const blob = await res.blob();
    const url  = URL.createObjectURL(blob);
    const a = document.createElement("a");
    a.href = url; a.download = "nanotubes.zip"; a.click();
    URL.revokeObjectURL(url);
    $("batch-overlay").classList.remove("open");
    toast(t("batch.done", { n: chiralities.length }), "ok", 5000);
  } catch (err) {
    toast(t("batch.fail") + ": " + err.message, "err", 6000);
  } finally {
    btn.disabled = false;
    btn.textContent = "";
    // Restaurar via i18n (suporta innerHTML para o &amp;)
    btn.innerHTML = t("batch.go");
  }
}

/* ---------- Init ----------------------------------------------------------- */
// Aplicar idioma persistido antes de iniciar lógica que dependa de t()
(function initLang() {
  let saved = "pt";
  try { saved = localStorage.getItem("ntbuilder-lang") || "pt"; } catch (_) {}
  applyLang(saved);
})();
loadExamples();
renderCutoffs();
updateActionButtons();
refreshHistoryButtons();
