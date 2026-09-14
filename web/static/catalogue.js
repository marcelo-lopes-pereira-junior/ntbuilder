/* Catálogo de nanotubos — a página.
 *
 * Quem conta é o servidor.  A contagem cai no cubo de agregados e responde em
 * milissegundos, então cada mexida refaz o número e a amostra (com um freio de
 * 250 ms).  Os controles deslizantes andam nas bordas das faixas do cubo, e é
 * isso que permite prometer um número exato.
 *
 * A unidade é o tubo distinto.  Um (n,m) de um Sistema 2D assimétrico vale
 * dois tubos, um por sentido de enrolamento, e o servidor já soma assim.  A
 * amostra mostra uma linha por (n,m), com os sentidos que passam no filtro.
 */
"use strict";

const API = "api/cat";   // relativa, porque o site vive atrás de /ntbuilder
const PAGE = 12;         // linhas por página da amostra
const INSTANT_MAX_SECONDS = 60;

/* ── Textos ─────────────────────────────────────────────────────────────── */
const I18N = {
  pt: {
    "nav.builder": "Construtor", "nav.catalogue": "Catálogo", "nav.tools": "Ferramentas",
    "hero.title": "Catálogo de nanotubos",
    "hero.lead": "Para cada Sistema 2D de oito Bases de Dados 2D públicas, sem estruturas repetidas entre elas, o catálogo reúne todos os nanotubos distintos que se pode enrolar a partir dele, com a menor célula unitária que mantém a periodicidade dentro de 0,5 % e a verificação de ligações espúrias. Escolha as bases ao lado, marque os elementos e ajuste o filtro, e a contagem e a amostra acompanham cada mudança. Como tudo é calculado está explicado no fim da página.",
    "hero.more": "Como o catálogo funciona",
    "db.title": "Base de Dados 2D", "db.systems": "Sistemas 2D",
    "stat.systems": "Sistemas 2D", "stat.tubes": "Nanotubos",
    "stat.exact": "Com célula exata", "stat.clean": "Sem ligação espúria",
    "stat.pct": "{p} % dos nanotubos",
    "pt.title": "Elementos", "pt.in": "Contém", "pt.out": "Não contém",
    "pt.excludeSwitch": "Excluir elementos",
    "pt.help": "Digite os símbolos separados por espaço ou clique na tabela. O clique vai para o campo destacado.",
    "pt.dim": "Não combina com a seleção", "pt.absent": "Não está no banco",
    "pt.materials": "Sistemas 2D", "pt.withSel": "Sistemas 2D com a seleção",
    "pt.clear": "Limpar elementos",
    "pt.classesIn": "Incluir classes", "pt.classesOut": "Excluir classes",
    "pt.classesHint": "Excluir uma classe tira os Sistemas 2D que têm qualquer elemento dela. Para ficar só com alguns elementos, exclua todos e clique na tabela nos que quer manter.",
    "pt.cAll": "Todos", "pt.outMany": "{n} elementos excluídos",
    "pt.cTM": "Metais de transição", "pt.cLan": "Lantanídeos", "pt.cAct": "Actinídeos",
    "pt.cRad": "Radioativos", "pt.cAlk": "Metais alcalinos", "pt.cAlkE": "Metais alcalino-terrosos",
    "pt.cHal": "Halogênios",
    "pt.warnUnknown": "Símbolo desconhecido, {s}.",
    "pt.warnAbsent": "{s} não aparece em nenhum Sistema 2D do banco.",
    "pt.warnDim": "{s} não combina com a seleção atual.",
    "pt.warnAll": "{s} está em todos os Sistemas 2D da seleção. Excluí-lo deixaria a seleção vazia.",
    "filter.title": "Filtro", "filter.lattices": "Classe de rede", "filter.latHelp": "O que é",
    "filter.thickness": "Espessura da camada", "filter.thickHint": "Do átomo mais baixo ao mais alto da camada.",
    "filter.diameter": "Diâmetro", "filter.cavity": "Cavidade interna mínima",
    "filter.cavityHint": "D − t, entre os centros dos átomos da parede interna.",
    "filter.atoms": "Átomos por célula", "filter.eps": "Resíduo de periodicidade",
    "filter.epsExact": "Exato", "filter.chirality": "Quiralidade",
    "filter.achiral": "Aquiral", "filter.chiral": "Quiral",
    "filter.chirHint": "Aquiral significa Ch paralelo ou perpendicular a um espelho ou deslizamento do Sistema 2D.",
    "filter.bonds": "Ligações espúrias", "filter.clean": "Sem ligação espúria", "filter.cleanAll": "Todos os tubos",
    "filter.reset": "Limpar filtros",
    "range.to": "a", "range.upTo": "Até",
    "count.tubes": "Tubos na seleção", "count.atoms": "Átomos", "count.size": "Arquivo compactado",
    "count.time": "Construção", "count.goto": "Ir para o pedido",
    "badge.exact": "Contagem exata", "badge.approx": "Estimativa aproximada",
    "src.cubo": "Cubo", "src.indice": "Índice",
    "sample.title": "Amostra da seleção",
    "sample.note": "Cada linha é um Sistema 2D da seleção, com quantos tubos ele tem nela. Abra a linha para ver os tubos. Quando os dois sentidos de enrolamento dão tubos diferentes, cada tubo traz os arquivos dos sentidos que passam no filtro.",
    "sample.infoSys": "Sistemas 2D {a} a {b} de {n}", "sample.tubesInfo": "Tubos {a} a {b} de {n}",
    "sample.lattice": "Classe de rede", "sample.thick": "Espessura (Å)", "sample.tubes": "Tubos",
    "sample.diam": "D (Å)", "sample.minCell": "Menor célula (átomos)",
    "sample.expand": "Ver os tubos deste Sistema 2D", "sample.openLayerTitle": "Abre a camada 2D deste Sistema 2D no construtor, numa nova aba",
    "sample.info": "Linhas {a} a {b} de {n}",
    "sample.page": "Página", "sample.of": "de",
    "sample.first": "Primeira", "sample.prev": "Anterior", "sample.next": "Próxima", "sample.last": "Última",
    "sample.noRows": "Nenhum tubo nesta seleção. Solte algum filtro.",
    "sample.system": "Sistema 2D", "sample.cavity": "Cavidade (Å)", "sample.atoms": "Átomos",
    "sample.chir": "Quiralidade", "sample.bonds": "Espúrias", "sample.files": "Arquivos",
    "sample.none": "Nenhuma", "sample.download": "Baixar", "sp.formed": "formada", "sp.broken": "rompida",
    "sample.open": "Abrir no NTBuilder", "sample.openTitle": "Abre a camada 2D deste Sistema 2D no construtor, numa nova aba, com este (n,m)",
    "sample.ccw": "Anti-horário", "sample.cw": "Horário",
    "sample.ccwTitle": "Visto da ponta do eixo para onde T aponta, Ch avança no sentido anti-horário",
    "sample.cwTitle": "Visto da ponta do eixo para onde T aponta, Ch avança no sentido horário",
    "req.title": "Pedido",
    "req.lead": "O pedido leva a seleção definida acima. Nos Sistemas 2D em que os dois sentidos de enrolamento dão tubos diferentes, os dois vão no arquivo. Escolha o formato, confira o resumo e confirme com nome e e-mail.",
    "req.formats": "Formato", "req.vacuum": "Vácuo (Å)", "req.perMat": "Máximo por Sistema 2D",
    "req.perMatHint": "Zero inclui todos. A contagem acima não aplica este limite.",
    "req.name": "Nome", "req.email": "E-mail", "req.send": "Confirmar pedido",
    "req.note": "Seleções com construção estimada em até 1 minuto são montadas na hora no servidor, e o download começa na tela assim que o arquivo fica pronto. As maiores vão para o cluster do laboratório, e o link chega por e-mail quando terminam. Nos dois casos o protocolo e o link, válido por 7 dias, vão para o seu e-mail.",
    "req.sumTubes": "Arquivos de tubo", "req.sumAtoms": "Átomos", "req.sumSize": "Arquivo compactado",
    "req.sumTime": "Construção", "req.sumFmt": "Formato", "req.sumMode": "Entrega",
    "req.modeNow": "Download na hora", "req.modeMail": "Link por e-mail",
    "track.title": "Acompanhar um pedido", "track.check": "Consultar", "track.placeholder": "Número do protocolo",
    "state.queued": "Na fila do cluster", "state.localQueued": "Aguardando a vez no servidor",
    "state.running": "Construindo", "state.done": "Pronto",
    "state.failed": "Falhou", "state.expired": "Expirado",
    "build.title": "Montando o arquivo", "build.wait": "Aguardando a vez no servidor",
    "build.running": "Construindo {d} de {n} tubos", "build.done": "Pronto. O download começou.",
    "build.failed": "A montagem falhou. {m}",
    "build.note": "O protocolo e o link, válido por 7 dias, também vão para o seu e-mail.",
    "build.again": "Baixar de novo", "build.close": "Fechar",
    "msg.empty": "A seleção está vazia. Solte algum filtro.",
    "msg.needEmail": "Preencha um e-mail válido.",
    "msg.sent": "Pedido registrado. O protocolo está na tela e no seu e-mail.",
    "msg.instant": "Pedido registrado. O arquivo está sendo montado.",
    "msg.mailFail": "Pedido registrado, mas o e-mail não saiu. Guarde o protocolo.",
    "msg.cached": "Esta seleção já estava pronta. O arquivo é o mesmo.",
    "msg.noProto": "Protocolo não encontrado.",
    "dl.link": "Baixar", "dl.expires": "Expira em",
    "fig.title": "Classes de rede",
    "fig.caption": "As cinco redes de Bravais bidimensionais e como a métrica de cada uma define os tubos que ela admite.",
    "ex.title": "Como o catálogo funciona",
    "ex.src.t": "De onde vêm os Sistemas 2D",
    "ex.src.p": "As estruturas relaxadas vêm de oito Bases de Dados 2D, o C2DB, o 2DMatPedia, as redes de carbono 2D do material suplementar aberto de um artigo do J. Phys. Chem. Lett. de 2021, o MC2D, o JARVIS-DFT 2D, a base de MXenes, o MatHub-2d e o Alexandria, este até 0,1 eV por átomo acima do casco convexo. Uma estrutura presente em mais de uma base entra uma vez só. Duas estruturas contam como a mesma quando têm a mesma fórmula, área por átomo e espessura próximas, o mesmo número de átomos na célula primitiva e, sobrepostas, nenhum átomo longe do seu par. Antes de enrolar, cada estrutura passa pela busca de simetria da ferramenta, com uma escada de tolerâncias, para que pequenas distorções numéricas não escondam a simetria real da rede.",
    "ex.nm.t": "Quais tubos entram",
    "ex.nm.p": "Um tubo é definido pelo vetor quiral Ch = n·a₁ + m·a₂. Pares (n,m) que a simetria do Sistema 2D leva uns nos outros dão o mesmo tubo, e por isso entra um representante de cada família. O diâmetro D = |Ch|/π fica entre a espessura t da camada e o maior valor entre 30 Å e t + 30 Å. Abaixo de t a parede atravessaria o eixo, e o teto acompanha a espessura para que camadas grossas também tenham tubos com cavidade.",
    "ex.d.t": "Diâmetro e cavidade",
    "ex.d.p": "D é medido no meio da espessura. Cada átomo fica no raio D/2 + z, com z contado a partir do plano médio da camada, e a parede interna fica em D/2 − t/2. A cavidade interna D − t é o diâmetro do círculo que passa pelos centros dos átomos da parede interna, sem descontar o raio atômico.",
    "ex.cell.t": "A célula de cada tubo",
    "ex.cell.p": "O vetor de translação T precisa ser um vetor da rede perpendicular a Ch. A célula exata existe quando os parâmetros da rede são comensuráveis, e o seu tamanho depende de quão simples é essa razão. Em redes hexagonais ela é sempre pequena, e em redes de razão quase irracional pode passar de milhões de átomos. Guardamos a menor célula cujo resíduo de periodicidade ε não passa de 0,5 %, com até 50 000 átomos. O resíduo equivale a um cisalhamento uniforme de ordem ε, e ε = 0 marca a célula exata.",
    "ex.sense.t": "Sentido de enrolamento",
    "ex.sense.p": "Uma camada pode ser enrolada com qualquer uma das duas faces para fora. Quando a camada virada coincide com a original por alguma operação de simetria, como no grafeno, no MoS₂ e no bifenileno, os dois sentidos dão o mesmo tubo e o catálogo guarda um só. Quando não coincide, como nos Janus e no penta-grafeno, os dois tubos são diferentes e o catálogo guarda os dois. Eles são chamados de anti-horário e horário. Visto da ponta do eixo para onde T aponta, no anti-horário Ch avança no sentido anti-horário e os átomos de z positivo da camada, como gravada na Base de Dados 2D, ficam na parede externa, e no horário ficam na parede interna.",
    "ex.bonds.t": "Ligações espúrias",
    "ex.bonds.p": "Enrolar aproxima alguns átomos e afasta outros. Cada par de átomos do tubo é comparado com o mesmo par na camada plana, com o corte de 1,2 vez a soma dos raios covalentes. Um par sem ligação no Sistema 2D que fica abaixo do corte e ao menos 10 % mais perto é uma ligação formada, e um par ligado que passa do corte e fica ao menos 10 % mais longe é uma ligação rompida. As duas contam como ligação espúria. Com o filtro ligado só entram os tubos sem nenhuma, e nos Sistemas 2D com dois sentidos cada tubo é julgado por conta própria.",
    "ex.chir.t": "Quiralidade",
    "ex.chir.p": "Um tubo é aquiral quando Ch é paralelo ou perpendicular a uma linha de espelho ou de deslizamento do Sistema 2D, e nesse caso ele coincide com a sua imagem especular. No grafeno isso dá os tubos zigzag e armchair, no bifenileno dá (n,0) e (0,m), e no penta-grafeno dá as diagonais (n,n) e (n,−n). Os demais são quirais.",
    "ex.req.t": "O pedido",
    "ex.req.p": "A contagem da seleção é exata e responde na hora. Ao confirmar o pedido com nome e e-mail, você recebe um protocolo. Se a construção estimada leva até 1 minuto, o arquivo é montado no próprio servidor e o download começa na tela assim que fica pronto. Se leva mais, os tubos são construídos no cluster do laboratório e o link chega por e-mail. Nos dois casos o link fica válido por 7 dias. O arquivo traz um CIF por tubo, ou o formato escolhido, o manifest.csv com a origem e os números de cada tubo e o CITATION.txt com as referências das bases usadas.",
    "footer.cite": "Se você usa o catálogo, por favor cite o NTBuilder e as bases de origem.",
    "lat.hexagonal": "Hexagonal", "lat.square": "Quadrada", "lat.rectangular": "Retangular",
    "lat.centred rectangular": "Retangular centrada", "lat.oblique": "Oblíqua",
    "kind.0": "Aquiral", "kind.1": "Quiral",
  },
  en: {
    "nav.builder": "Builder", "nav.catalogue": "Catalogue", "nav.tools": "Tools",
    "hero.title": "Nanotube catalogue",
    "hero.lead": "For every 2D System of eight public 2D Databases, with no structure repeated across them, the catalogue gathers all distinct nanotubes that can be rolled from it, with the smallest unit cell that keeps periodicity within 0.5 % and the spurious-bond check. Choose the databases on the side, mark the elements and adjust the filter, and the count and the sample follow every change. How everything is computed is explained at the end of the page.",
    "hero.more": "How the catalogue works",
    "db.title": "2D Database", "db.systems": "2D Systems",
    "stat.systems": "2D Systems", "stat.tubes": "Nanotubes",
    "stat.exact": "With an exact cell", "stat.clean": "No spurious bond",
    "stat.pct": "{p} % of the nanotubes",
    "pt.title": "Elements", "pt.in": "Contains", "pt.out": "Excludes",
    "pt.excludeSwitch": "Exclude elements",
    "pt.help": "Type the symbols separated by spaces or click on the table. A click goes to the highlighted field.",
    "pt.dim": "Does not combine with the selection", "pt.absent": "Not in the database",
    "pt.materials": "2D Systems", "pt.withSel": "2D Systems with the selection",
    "pt.clear": "Clear elements",
    "pt.classesIn": "Include classes", "pt.classesOut": "Exclude classes",
    "pt.classesHint": "Excluding a class removes the 2D Systems that have any of its elements. To keep only some elements, exclude all and click on the table the ones you want to keep.",
    "pt.cAll": "All", "pt.outMany": "{n} elements excluded",
    "pt.cTM": "Transition metals", "pt.cLan": "Lanthanides", "pt.cAct": "Actinides",
    "pt.cRad": "Radioactive", "pt.cAlk": "Alkali metals", "pt.cAlkE": "Alkaline earth metals",
    "pt.cHal": "Halogens",
    "pt.warnUnknown": "Unknown symbol, {s}.",
    "pt.warnAbsent": "{s} does not appear in any 2D System of the database.",
    "pt.warnDim": "{s} does not combine with the current selection.",
    "pt.warnAll": "{s} is in every 2D System of the selection. Excluding it would empty the selection.",
    "filter.title": "Filter", "filter.lattices": "Lattice class", "filter.latHelp": "What is it",
    "filter.thickness": "Layer thickness", "filter.thickHint": "From the lowest to the highest atom of the layer.",
    "filter.diameter": "Diameter", "filter.cavity": "Minimum inner cavity",
    "filter.cavityHint": "D − t, between the centres of the inner-wall atoms.",
    "filter.atoms": "Atoms per cell", "filter.eps": "Periodicity residual",
    "filter.epsExact": "Exact", "filter.chirality": "Chirality",
    "filter.achiral": "Achiral", "filter.chiral": "Chiral",
    "filter.chirHint": "Achiral means Ch parallel or perpendicular to a mirror or glide of the 2D System.",
    "filter.bonds": "Spurious bonds", "filter.clean": "No spurious bond", "filter.cleanAll": "All tubes",
    "filter.reset": "Clear filters",
    "range.to": "to", "range.upTo": "Up to",
    "count.tubes": "Tubes in the selection", "count.atoms": "Atoms", "count.size": "Compressed file",
    "count.time": "Build time", "count.goto": "Go to the request",
    "badge.exact": "Exact count", "badge.approx": "Approximate estimate",
    "src.cubo": "Cube", "src.indice": "Index",
    "sample.title": "Sample of the selection",
    "sample.note": "Each row is one 2D System of the selection, with how many of its tubes are in it. Open the row to see the tubes. When the two rolling senses give different tubes, each tube carries the files of the senses that pass the filter.",
    "sample.infoSys": "2D Systems {a} to {b} of {n}", "sample.tubesInfo": "Tubes {a} to {b} of {n}",
    "sample.lattice": "Lattice class", "sample.thick": "Thickness (Å)", "sample.tubes": "Tubes",
    "sample.diam": "D (Å)", "sample.minCell": "Smallest cell (atoms)",
    "sample.expand": "See the tubes of this 2D System", "sample.openLayerTitle": "Opens the 2D layer of this 2D System in the builder, in a new tab",
    "sample.info": "Rows {a} to {b} of {n}",
    "sample.page": "Page", "sample.of": "of",
    "sample.first": "First", "sample.prev": "Previous", "sample.next": "Next", "sample.last": "Last",
    "sample.noRows": "No tube in this selection. Loosen a filter.",
    "sample.system": "2D System", "sample.cavity": "Cavity (Å)", "sample.atoms": "Atoms",
    "sample.chir": "Chirality", "sample.bonds": "Spurious", "sample.files": "Files",
    "sample.none": "None", "sample.download": "Download", "sp.formed": "formed", "sp.broken": "broken",
    "sample.open": "Open in NTBuilder", "sample.openTitle": "Opens the 2D layer of this 2D System in the builder, in a new tab, with this (n,m)",
    "sample.ccw": "Counterclockwise", "sample.cw": "Clockwise",
    "sample.ccwTitle": "Seen from the end of the axis T points to, Ch advances counterclockwise",
    "sample.cwTitle": "Seen from the end of the axis T points to, Ch advances clockwise",
    "req.title": "Request",
    "req.lead": "The request takes the selection defined above. For 2D Systems whose two rolling senses give different tubes, both go into the file. Choose the format, check the summary and confirm with name and e-mail.",
    "req.formats": "Format", "req.vacuum": "Vacuum (Å)", "req.perMat": "Maximum per 2D System",
    "req.perMatHint": "Zero includes all. The count above does not apply this limit.",
    "req.name": "Name", "req.email": "E-mail", "req.send": "Confirm request",
    "req.note": "Selections with an estimated build time of up to 1 minute are assembled right away on the server, and the download starts on screen as soon as the file is ready. Larger ones go to the laboratory cluster, and the link arrives by e-mail when they finish. In both cases the protocol and the link, valid for 7 days, are sent to your e-mail.",
    "req.sumTubes": "Tube files", "req.sumAtoms": "Atoms", "req.sumSize": "Compressed file",
    "req.sumTime": "Build time", "req.sumFmt": "Format", "req.sumMode": "Delivery",
    "req.modeNow": "Download right away", "req.modeMail": "Link by e-mail",
    "track.title": "Track a request", "track.check": "Check", "track.placeholder": "Protocol number",
    "state.queued": "Queued on the cluster", "state.localQueued": "Waiting for its turn on the server",
    "state.running": "Building", "state.done": "Ready",
    "state.failed": "Failed", "state.expired": "Expired",
    "build.title": "Assembling the file", "build.wait": "Waiting for its turn on the server",
    "build.running": "Building {d} of {n} tubes", "build.done": "Ready. The download has started.",
    "build.failed": "The build failed. {m}",
    "build.note": "The protocol and the link, valid for 7 days, are also sent to your e-mail.",
    "build.again": "Download again", "build.close": "Close",
    "msg.empty": "The selection is empty. Loosen a filter.",
    "msg.needEmail": "Please give a valid e-mail.",
    "msg.sent": "Request registered. The protocol is on screen and in your e-mail.",
    "msg.instant": "Request registered. The file is being assembled.",
    "msg.mailFail": "Request registered, but the e-mail did not go out. Keep the protocol.",
    "msg.cached": "This selection was already built. Same file.",
    "msg.noProto": "Protocol not found.",
    "dl.link": "Download", "dl.expires": "Expires on",
    "fig.title": "Lattice classes",
    "fig.caption": "The five two-dimensional Bravais lattices and how the metric of each defines the tubes it admits.",
    "ex.title": "How the catalogue works",
    "ex.src.t": "Where the 2D Systems come from",
    "ex.src.p": "The relaxed structures come from eight 2D Databases, C2DB, 2DMatPedia, the 2D carbon networks in the open supplementary material of a 2021 J. Phys. Chem. Lett. article, MC2D, JARVIS-DFT 2D, the MXene database, MatHub-2d and Alexandria, the latter up to 0.1 eV per atom above the convex hull. A structure present in more than one database enters once. Two structures count as the same when they share the formula, have close area per atom and thickness, the same number of atoms in the primitive cell and, once overlaid, no atom far from its partner. Before rolling, every structure goes through the tool's symmetry search, with a ladder of tolerances, so that small numerical distortions do not hide the true symmetry of the lattice.",
    "ex.nm.t": "Which tubes are included",
    "ex.nm.p": "A tube is defined by the chiral vector Ch = n·a₁ + m·a₂. Pairs (n,m) that the symmetry of the 2D System maps onto each other give the same tube, so one representative of each family is kept. The diameter D = |Ch|/π lies between the layer thickness t and the larger of 30 Å and t + 30 Å. Below t the wall would cross the axis, and the ceiling follows the thickness so that thick layers also have tubes with a cavity.",
    "ex.d.t": "Diameter and cavity",
    "ex.d.p": "D is measured at the middle of the thickness. Each atom sits at radius D/2 + z, with z counted from the mid-plane of the layer, and the inner wall lies at D/2 − t/2. The inner cavity D − t is the diameter of the circle through the centres of the inner-wall atoms, without subtracting atomic radii.",
    "ex.cell.t": "The cell of each tube",
    "ex.cell.p": "The translation vector T must be a lattice vector perpendicular to Ch. The exact cell exists when the lattice parameters are commensurate, and its size depends on how simple that ratio is. On hexagonal lattices it is always small, and on lattices with a nearly irrational ratio it can exceed millions of atoms. We keep the smallest cell whose periodicity residual ε does not exceed 0.5 %, with up to 50,000 atoms. The residual amounts to a uniform shear of order ε, and ε = 0 marks the exact cell.",
    "ex.sense.t": "Rolling sense",
    "ex.sense.p": "A layer can be rolled with either face outside. When the flipped layer coincides with the original through some symmetry operation, as in graphene, MoS₂ and biphenylene, both senses give the same tube and the catalogue keeps one. When it does not, as in Janus layers and penta-graphene, the two tubes are different and the catalogue keeps both. They are called counterclockwise and clockwise. Seen from the end of the axis T points to, in the counterclockwise tube Ch advances counterclockwise and the atoms with positive z in the layer, as stored in the 2D Database, sit on the outer wall, and in the clockwise tube they sit on the inner wall.",
    "ex.bonds.t": "Spurious bonds",
    "ex.bonds.p": "Rolling brings some atoms together and pulls others apart. Every pair of atoms of the tube is compared with the same pair in the flat layer, with a cutoff of 1.2 times the sum of covalent radii. A pair not bonded in the 2D System that falls below the cutoff and gets at least 10 % closer is a formed bond, and a bonded pair that passes the cutoff and gets at least 10 % farther is a broken bond. Both count as spurious bonds. With the filter on only tubes without any are included, and for 2D Systems with two senses each tube is judged on its own.",
    "ex.chir.t": "Chirality",
    "ex.chir.p": "A tube is achiral when Ch is parallel or perpendicular to a mirror or glide line of the 2D System, and in that case it coincides with its mirror image. In graphene this gives the zigzag and armchair tubes, in biphenylene (n,0) and (0,m), and in penta-graphene the diagonals (n,n) and (n,−n). All others are chiral.",
    "ex.req.t": "The request",
    "ex.req.p": "The selection count is exact and answers immediately. When you confirm the request with name and e-mail, you get a protocol. If the estimated build takes up to 1 minute, the file is assembled on the server itself and the download starts on screen as soon as it is ready. If it takes longer, the tubes are built on the laboratory cluster and the link arrives by e-mail. In both cases the link stays valid for 7 days. The file has one CIF per tube, or the chosen format, manifest.csv with the origin and numbers of every tube, and CITATION.txt with the references of the databases used.",
    "footer.cite": "If you use the catalogue, please cite NTBuilder and the source databases.",
    "lat.hexagonal": "Hexagonal", "lat.square": "Square", "lat.rectangular": "Rectangular",
    "lat.centred rectangular": "Centred rectangular", "lat.oblique": "Oblique",
    "kind.0": "Achiral", "kind.1": "Chiral",
  },
};
let LANG = "pt";
const t = (k, vars) => {
  let s = (I18N[LANG] && I18N[LANG][k]) || I18N.pt[k] || k;
  if (vars) Object.entries(vars).forEach(([a, b]) => { s = s.replace(`{${a}}`, b); });
  return s;
};

const SOURCE_LABEL = { c2db: "C2DB", "2dmatpedia": "2DMatPedia", jz1c03193: "Redes de carbono 2D", mc2d: "MC2D", jarvis: "JARVIS-DFT 2D", alexandria: "Alexandria 2D", mxene: "MXenes", mathub2d: "MatHub-2d" };
const SOURCE_LABEL_EN = { jz1c03193: "2D carbon networks" };
const sourceLabel = (s) => (LANG === "en" && SOURCE_LABEL_EN[s]) || SOURCE_LABEL[s] || s;

/* ── Formatação ─────────────────────────────────────────────────────────── */
const locale = () => (LANG === "pt" ? "pt-BR" : "en-US");
const fmtInt = (n) => (n == null ? "—" : Number(n).toLocaleString(locale()));
const fmtNum = (x, d) => Number(x).toLocaleString(locale(), { minimumFractionDigits: d, maximumFractionDigits: d });
function fmtBytes(n) {
  if (n == null) return "—";
  const u = ["B", "kB", "MB", "GB", "TB"];
  let i = 0;
  while (n >= 1024 && i < u.length - 1) { n /= 1024; i += 1; }
  return `${fmtNum(n, i ? 1 : 0)} ${u[i]}`;
}
function fmtTime(s) {
  if (s == null) return "—";
  if (s < 90) return `${Math.max(1, Math.round(s))} s`;
  if (s < 5400) return `${Math.round(s / 60)} min`;
  return `${fmtNum(s / 3600, 1)} h`;
}
function toast(msg, kind) {
  const box = document.getElementById("toast-container");
  if (!box) return;
  const el = document.createElement("div");
  el.className = `toast${kind ? ` ${kind}` : ""}`;
  el.textContent = msg;
  box.appendChild(el);
  setTimeout(() => el.remove(), 6000);
}
async function post(path, body, query) {
  const r = await fetch(`${API}${path}${query || ""}`, {
    method: "POST", headers: { "Content-Type": "application/json" }, body: JSON.stringify(body),
  });
  const data = await r.json().catch(() => ({}));
  if (!r.ok) throw new Error(typeof data.detail === "string" ? data.detail : `Erro ${r.status}`);
  return data;
}

/* ── Estado ─────────────────────────────────────────────────────────────── */
const state = {
  meta: null, count: null,
  elementCounts: {},     // Sistemas 2D por elemento, no banco todo
  liveCounts: null,      // o mesmo com a seleção atual (null = sem restrição)
  liveTotal: null,       // quantos Sistemas 2D a seleção deixa
  mode: "in",            // campo que recebe o clique na tabela
  classIn: new Set(),    // classes incluídas (ao menos um elemento de cada)
  dMax: 80, trackTimer: null,
  page: 1, sampleRows: [], sampleFilter: null,
  grouped: false, systems: [], expanded: new Map(),   // amostra por Sistema 2D
  build: null, buildTimer: null,
};
const elState = new Map();   // símbolo -> "in" | "out"

const T_BINS = [0, 1, 2, 3, 4, 5, 6, 8, 10, 15, 20, 30, 60];
const DIN_STEPS = [0, 5, 10, 20];
const ATOM_STEPS = [100, 200, 500, 1000, 2000, 5000, 10000, 20000, 50000];
const EPS_STEPS = ["exact", 0.1, 0.2, 0.5];
const $ = (id) => document.getElementById(id);

/* ── Controles ──────────────────────────────────────────────────────────── */
const onValues = (id) => Array.from(document.querySelectorAll(`#${id} .is-on`)).map((b) => b.dataset.v);

function wirePills(id, onChange) {
  $(id).addEventListener("click", (ev) => {
    const b = ev.target.closest("[data-v]");
    if (!b) return;
    b.classList.toggle("is-on");
    b.setAttribute("aria-pressed", String(b.classList.contains("is-on")));
    onChange();
  });
}
function wireSeg(id, onChange) {
  $(id).addEventListener("click", (ev) => {
    const b = ev.target.closest("[data-v]");
    if (!b) return;
    document.querySelectorAll(`#${id} [data-v]`).forEach((x) => x.classList.toggle("is-on", x === b));
    onChange();
  });
}

function fillBetween(fillId, lo, hi, max) {
  const f = $(fillId);
  f.style.left = `${(lo / max) * 100}%`;
  f.style.right = `${100 - (hi / max) * 100}%`;
}

function sliderReadouts() {
  const to = t("range.to");
  const dmin = +$("f-d-min").value, dmax = +$("f-d-max").value;
  $("v-diam").textContent = `${dmin} ${to} ${dmax} Å`;
  fillBetween("f-diam-fill", dmin, dmax, +$("f-d-max").max);
  const tlo = +$("f-t-min").value, thi = +$("f-t-max").value;
  $("v-thick").textContent = `${fmtInt(T_BINS[tlo])} ${to} ${fmtInt(T_BINS[thi])} Å`;
  fillBetween("f-thick-fill", tlo, thi, T_BINS.length - 1);
  $("v-din").textContent = `${DIN_STEPS[+$("f-din").value]} Å`;
  $("v-atoms").textContent = `${t("range.upTo")} ${fmtInt(ATOM_STEPS[+$("f-atoms").value])}`;
  const e = EPS_STEPS[+$("f-eps").value];
  $("v-eps").textContent = e === "exact" ? t("filter.epsExact") : `${t("range.upTo")} ${fmtNum(e, 1)} %`;
}

function wireDual(loId, hiId, onChange) {
  const lo = $(loId), hi = $(hiId);
  lo.addEventListener("input", () => { if (+lo.value > +hi.value) lo.value = hi.value; sliderReadouts(); onChange(); });
  hi.addEventListener("input", () => { if (+hi.value < +lo.value) hi.value = lo.value; sliderReadouts(); onChange(); });
}

/** O filtro como o servidor o espera: uma função só, para a contagem, a
 *  amostra e o pedido nunca perguntarem coisas diferentes. */
function readFilter() {
  const all = (got, total) => (total && got.length >= total ? null : got);
  const eps = EPS_STEPS[+$("f-eps").value];
  const chir = onValues("f-chirality");
  const list = (kind) => { const v = [...elState].filter(([, k]) => k === kind).map(([s]) => s); return v.length ? v : null; };
  const groups = CLASSES.filter((c) => state.classIn.has(c.key)).map((c) => classPresent(c)).filter((g) => g.length);
  return {
    sources: all(onValues("f-sources"), state.meta ? state.meta.by_source.length : 0),
    lattices: all(onValues("f-lattices"), state.meta ? state.meta.by_lattice.length : 0),
    elements_all: list("in"),
    elements_any: groups.length ? groups : null,
    elements_none: list("out"),
    t_min: T_BINS[+$("f-t-min").value],
    t_max: T_BINS[+$("f-t-max").value],
    d_min: +$("f-d-min").value,
    d_max: +$("f-d-max").value,
    din_min: DIN_STEPS[+$("f-din").value],
    atoms_max: ATOM_STEPS[+$("f-atoms").value],
    eps_max: eps === "exact" ? 0.5 : eps,
    exact_only: eps === "exact",
    chirality: chir.length === 2 ? null : chir,
    clean: onValues("f-clean")[0] === "all" ? "all" : "clean",
    per_material_max: Math.max(0, parseInt($("f-per-mat").value, 10) || 0),
    formats: onValues("f-formats").length ? onValues("f-formats") : ["cif"],
    vacuum: +$("f-vacuum").value || 10,
  };
}
function filterIsEmpty() {
  return ["f-sources", "f-lattices", "f-chirality", "f-formats"].some((id) => onValues(id).length === 0);
}

/* ── Números por Base de Dados 2D ───────────────────────────────────────── */
function renderStats() {
  const on = new Set(onValues("f-sources"));
  const sum = { n_materials: 0, n_tubes: 0, n_exact: 0, n_clean: 0 };
  (state.meta.by_source_stats || []).forEach((s) => {
    if (on.has(s.source)) Object.keys(sum).forEach((k) => { sum[k] += s[k] || 0; });
  });
  const pct = (x) => (sum.n_tubes ? t("stat.pct", { p: fmtNum((100 * x) / sum.n_tubes, 1) }) : "");
  $("stat-systems").textContent = fmtInt(sum.n_materials);
  $("stat-tubes").textContent = fmtInt(sum.n_tubes);
  $("stat-exact").textContent = fmtInt(sum.n_exact);
  $("stat-exact-p").textContent = pct(sum.n_exact);
  $("stat-clean").textContent = fmtInt(sum.n_clean);
  $("stat-clean-p").textContent = pct(sum.n_clean);
  document.querySelectorAll("#f-sources .db-toggle").forEach((b) => {
    const s = (state.meta.by_source_stats || []).find((x) => x.source === b.dataset.v);
    b.querySelector(".db-toggle__name").textContent = sourceLabel(b.dataset.v);
    b.querySelector(".db-toggle__n").textContent = s ? `${fmtInt(s.n_materials)} ${t("db.systems")}` : "";
  });
}

/* ── Contagem ───────────────────────────────────────────────────────────── */
let countTimer = null, countSeq = 0;
function scheduleCount() { clearTimeout(countTimer); countTimer = setTimeout(refreshCount, 250); }

async function refreshCount() {
  refreshElementCounts();
  loadPage(1);
  if (filterIsEmpty()) {
    state.count = { n_tubes: 0, n_atoms: 0, bytes_zip: 0, seconds: 0, exact_count: true, n_rows: 0 };
    renderCount();
    return;
  }
  const seq = (countSeq += 1);
  $("cat-count").classList.add("is-busy");
  try {
    const data = await post("/count", readFilter(), "?fast=true");
    if (seq !== countSeq) return;
    state.count = data;
  } catch (e) {
    if (seq !== countSeq) return;
    state.count = null;
    toast(e.message, "err");
  } finally {
    if (seq === countSeq) { $("cat-count").classList.remove("is-busy"); renderCount(); }
  }
}

function renderCount() {
  const c = state.count;
  const set = (id, v) => { $(id).textContent = v; };
  if (!c) {
    ["c-tubes", "c-atoms", "c-size", "c-time"].forEach((i) => set(i, "—"));
    renderSummary(); renderPager(); return;
  }
  set("c-tubes", fmtInt(c.n_tubes));
  set("c-atoms", fmtInt(c.n_atoms));
  set("c-size", `~${fmtBytes(c.bytes_zip)}`);
  set("c-time", `~${fmtTime(c.seconds)}`);
  const badge = $("c-exact-badge");
  badge.textContent = c.exact_count ? t("badge.exact") : t("badge.approx");
  badge.className = `cat-badge${c.exact_count ? "" : " cat-badge--warn"}`;
  set("c-ms", c.ms != null ? `${fmtNum(c.ms, 0)} ms · ${t(`src.${c.source || "indice"}`)}` : "");
  $("req-send").disabled = !c.n_tubes;
  renderSummary();
  renderPager();
}

function renderSummary() {
  const c = state.count;
  const f = readFilter();
  const tile = (n, l) => `<div><span class="cat-tile__n">${n}</span><span class="cat-tile__l">${l}</span></div>`;
  $("req-summary").innerHTML = c
    ? tile(fmtInt((c.n_tubes || 0) * f.formats.length), t("req.sumTubes")) +
      tile(fmtInt(c.n_atoms), t("req.sumAtoms")) +
      tile(`~${fmtBytes(c.bytes_zip)}`, t("req.sumSize")) +
      tile(`~${fmtTime(c.seconds)}`, t("req.sumTime")) +
      tile(f.formats.map((x) => x.toUpperCase()).join(", "), t("req.sumFmt")) +
      tile((c.seconds || 0) <= (state.instantMax || INSTANT_MAX_SECONDS) ? t("req.modeNow") : t("req.modeMail"), t("req.sumMode"))
    : "";
}

/* ── Amostra paginada ───────────────────────────────────────────────────── */
let sampleSeq = 0;
function totalPages() {
  const n = state.grouped ? state.liveTotal : (state.count && state.count.n_rows);
  return n == null ? null : Math.max(1, Math.ceil(n / PAGE));
}

async function loadSample(page) {
  const seq = (sampleSeq += 1);
  if (filterIsEmpty()) {
    state.sampleRows = []; state.page = 1; state.sampleFilter = readFilter();
    renderSample();
    return;
  }
  const f = readFilter();
  const total = totalPages();
  page = Math.max(1, Math.floor(page) || 1);
  if (total && page > total && seq > 1) page = total;
  $("cat-sample").classList.add("is-busy");
  try {
    const data = await post("/sample", f, `?limit=${PAGE}&offset=${(page - 1) * PAGE}`);
    if (seq !== sampleSeq) return;
    state.page = page;
    state.sampleRows = data.tubes || [];
    state.sampleFilter = f;
  } catch (e) {
    if (seq !== sampleSeq) return;
    state.sampleRows = [];
    toast(e.message, "err");
  } finally {
    if (seq === sampleSeq) { $("cat-sample").classList.remove("is-busy"); renderSample(); }
  }
}

function spText(sp) {
  // "+A-B" ligação que o enrolamento forma, "-A-B" ligação que ele rompe.
  if (!sp) return t("sample.none");
  return sp.split(",").filter(Boolean).map((x) =>
    x[0] === "+" ? `${x.slice(1)} ${t("sp.formed")}` : x[0] === "-" ? `${x.slice(1)} ${t("sp.broken")}` : x).join(", ");
}

function sensesOf(x, cleanOnly) {
  // Um sentido só quando os dois dão o mesmo tubo, e aí a linha só entra no
  // filtro limpo se as duas checagens concordam.
  if (x.senses === 1) {
    return [{ key: null, inward: false, sp: x.spurious_ccw || x.spurious_cw }];
  }
  const both = [
    { key: "ccw", inward: false, clean: x.clean_ccw, sp: x.spurious_ccw },
    { key: "cw", inward: true, clean: x.clean_cw, sp: x.spurious_cw },
  ];
  return cleanOnly ? both.filter((s) => s.clean) : both;
}

function renderSample() {
  const f = state.sampleFilter || readFilter();
  const cleanOnly = f.clean === "clean";
  $("sample-table").classList.toggle("is-clean", cleanOnly);
  const body = document.querySelector("#sample-table tbody");
  body.innerHTML = "";
  if (!state.sampleRows.length) {
    body.innerHTML = `<tr><td colspan="9" class="cat-empty">${t("sample.noRows")}</td></tr>`;
    renderPager();
    return;
  }
  const fmt = f.formats[0] || "cif";
  for (const x of state.sampleRows) {
    const shown = sensesOf(x, cleanOnly);
    const link = (s) => `<a class="cat-link" ${s.key ? `title="${t(`sample.${s.key}Title`)}"` : ""} href="${API}/tube/${x.material_id}/${x.n}/${x.m}?fmt=${fmt}&vacuum=${f.vacuum}&inward=${s.inward}" download>${s.key ? t(`sample.${s.key}`) : t("sample.download")}</a>`;
    const sp = shown.map((s) => `${s.key ? `${t(`sample.${s.key}`)} ` : ""}${spText(s.sp)}`).join(" · ");
    const tr = document.createElement("tr");
    tr.innerHTML = `
      <td>${x.formula}<span class="cat-uid">${sourceLabel(x.source)} · ${x.uid}</span>${state.canOpen
        ? `<a class="cat-link cat-open" href="./?catalogo=${x.material_id}&n=${x.n}&m=${x.m}" target="_blank" rel="noopener" title="${t("sample.openTitle")}">${t("sample.open")} ↗</a>` : ""}</td>
      <td>(${x.n},${x.m})</td>
      <td>${fmtNum(x.diameter, 2)}</td>
      <td>${fmtNum(x.d_inner, 2)}</td>
      <td>${fmtInt(x.atoms)}</td>
      <td>${x.exact ? "0" : fmtNum(x.eps, 3)}</td>
      <td>${t(`kind.${x.kind}`)}</td>
      <td class="col-sp cat-sp">${sp}</td>
      <td>${shown.map(link).join(" · ")}</td>`;
    body.appendChild(tr);
  }
  renderPager();
}

function renderPager() {
  const total = totalPages();
  const p = state.page;
  const rows = state.grouped ? state.systems.length : state.sampleRows.length;
  $("pg-input").value = p;
  if (total) $("pg-input").max = total; else $("pg-input").removeAttribute("max");
  $("pg-total").textContent = total ? fmtInt(total) : "…";
  $("pg-first").disabled = p <= 1;
  $("pg-prev").disabled = p <= 1;
  const atEnd = total ? p >= total : rows < PAGE;
  $("pg-next").disabled = atEnd;
  $("pg-last").disabled = atEnd || !total;
  const n = state.grouped ? state.liveTotal : (state.count && state.count.n_rows);
  const a = rows ? (p - 1) * PAGE + 1 : 0;
  $("sample-info").textContent = n != null && rows
    ? t(state.grouped ? "sample.infoSys" : "sample.info", { a: fmtInt(a), b: fmtInt(a + rows - 1), n: fmtInt(n) }) : "";
}

function wirePager() {
  $("pg-first").addEventListener("click", () => loadPage(1));
  $("pg-prev").addEventListener("click", () => loadPage(state.page - 1));
  $("pg-next").addEventListener("click", () => loadPage(state.page + 1));
  $("pg-last").addEventListener("click", () => { const n = totalPages(); if (n) loadPage(n); });
  $("pg-input").addEventListener("change", () => {
    const n = totalPages();
    let p = parseInt($("pg-input").value, 10) || 1;
    if (n) p = Math.min(Math.max(1, p), n);
    loadPage(p);
  });
}

/* ── Amostra agrupada por Sistema 2D ────────────────────────────────────────
 * Uma linha por Sistema 2D, com os números dos seus tubos na seleção; um
 * clique abre os tubos dele, com páginas próprias.  Com Li marcado, a amostra
 * por tubo mostrava 35 páginas seguidas de Li2VF6. */
const SUB = 10;   // tubos por página dentro de um Sistema 2D

const loadPage = (page) => (state.grouped ? loadSystems(page) : loadSample(page));
const renderCurrent = () => (state.grouped ? renderSystems() : renderSample());

let sysSeq = 0;
async function loadSystems(page) {
  const seq = (sysSeq += 1);
  const f = readFilter();
  if (filterIsEmpty()) {
    state.systems = []; state.page = 1; state.sampleFilter = f; state.expanded.clear();
    renderSystems();
    return;
  }
  const total = totalPages();
  page = Math.max(1, Math.floor(page) || 1);
  if (total && page > total) page = total;
  $("cat-sample").classList.add("is-busy");
  try {
    const data = await post("/systems", f, `?limit=${PAGE}&offset=${(page - 1) * PAGE}`);
    if (seq !== sysSeq) return;
    state.page = page;
    state.systems = data.systems || [];
    state.sampleFilter = f;
    state.expanded.clear();
  } catch (e) {
    if (seq !== sysSeq) return;
    state.systems = [];
    toast(e.message, "err");
  } finally {
    if (seq === sysSeq) { $("cat-sample").classList.remove("is-busy"); renderSystems(); }
  }
}

function tubesHead() {
  return `<tr><th>(n,m)</th><th>D (Å)</th><th>${t("sample.cavity")}</th><th>${t("sample.atoms")}</th>`
    + `<th>ε (%)</th><th>${t("sample.chir")}</th><th class="col-sp">${t("sample.bonds")}</th><th>${t("sample.files")}</th></tr>`;
}

function renderSystems() {
  const f = state.sampleFilter || readFilter();
  const cleanOnly = f.clean === "clean";
  const table = $("sample-table");
  table.classList.toggle("is-clean", cleanOnly);
  table.classList.add("is-grouped");
  table.querySelector("thead").innerHTML = `<tr><th></th><th>${t("sample.system")}</th><th>${t("sample.lattice")}</th>`
    + `<th>${t("sample.thick")}</th><th>${t("sample.tubes")}</th><th>${t("sample.diam")}</th><th>${t("sample.minCell")}</th></tr>`;
  const body = table.querySelector("tbody");
  body.innerHTML = "";
  if (!state.systems.length) {
    body.innerHTML = `<tr><td colspan="7" class="cat-empty">${t("sample.noRows")}</td></tr>`;
    renderPager();
    return;
  }
  for (const x of state.systems) {
    const open = state.expanded.has(x.material_id);
    const tr = document.createElement("tr");
    tr.className = `cat-sys${open ? " is-open" : ""}`;
    tr.dataset.mid = x.material_id;
    tr.innerHTML = `
      <td><button type="button" class="cat-toggle" data-mid="${x.material_id}" aria-expanded="${open}" title="${t("sample.expand")}">${open ? "▾" : "▸"}</button></td>
      <td>${x.formula}<span class="cat-uid">${sourceLabel(x.source)} · ${x.uid}</span>${state.canOpen
        ? `<a class="cat-link cat-open" href="./?catalogo=${x.material_id}" target="_blank" rel="noopener" title="${t("sample.openLayerTitle")}">${t("sample.open")} ↗</a>` : ""}</td>
      <td>${t(`lat.${x.lattice_type}`)}</td>
      <td>${fmtNum(x.thickness, 2)}</td>
      <td>${fmtInt(x.n_tubes)}</td>
      <td>${fmtNum(x.d_min, 1)} ${t("range.to")} ${fmtNum(x.d_max, 1)}</td>
      <td>${fmtInt(x.atoms_min)}</td>`;
    body.appendChild(tr);
    if (open) {
      const sub = document.createElement("tr");
      sub.className = "cat-sys-tubes";
      sub.innerHTML = `<td colspan="7">${systemTubes(x, f, cleanOnly)}</td>`;
      body.appendChild(sub);
    }
  }
  renderPager();
}

function systemTubes(x, f, cleanOnly) {
  const ex = state.expanded.get(x.material_id);
  if (!ex || ex.loading) return `<div class="cat-sub cat-hint">…</div>`;
  const fmt = f.formats[0] || "cif";
  const rows = ex.rows.map((r) => {
    const shown = sensesOf(r, cleanOnly);
    const link = (s) => `<a class="cat-link" ${s.key ? `title="${t(`sample.${s.key}Title`)}"` : ""} href="${API}/tube/${r.material_id}/${r.n}/${r.m}?fmt=${fmt}&vacuum=${f.vacuum}&inward=${s.inward}" download>${s.key ? t(`sample.${s.key}`) : t("sample.download")}</a>`;
    const sp = shown.map((s) => `${s.key ? `${t(`sample.${s.key}`)} ` : ""}${spText(s.sp)}`).join(" · ");
    const open = state.canOpen ? ` · <a class="cat-link" href="./?catalogo=${r.material_id}&n=${r.n}&m=${r.m}" target="_blank" rel="noopener" title="${t("sample.openTitle")}">NTBuilder ↗</a>` : "";
    return `<tr><td>(${r.n},${r.m})</td><td>${fmtNum(r.diameter, 2)}</td><td>${fmtNum(r.d_inner, 2)}</td><td>${fmtInt(r.atoms)}</td>`
      + `<td>${r.exact ? "0" : fmtNum(r.eps, 3)}</td><td>${t(`kind.${r.kind}`)}</td><td class="col-sp cat-sp">${sp}</td>`
      + `<td>${shown.map(link).join(" · ")}${open}</td></tr>`;
  }).join("");
  const a = ex.rows.length ? (ex.page - 1) * SUB + 1 : 0;
  const b = (ex.page - 1) * SUB + ex.rows.length;
  const pages = Math.max(1, Math.ceil(x.n_rows / SUB));
  return `<div class="cat-sub">
      <table class="cat-table cat-table--sub${cleanOnly ? " is-clean" : ""}"><thead>${tubesHead()}</thead><tbody>${rows}</tbody></table>
      <div class="cat-sub__pager">
        <button type="button" data-sub="${x.material_id}" data-page="${ex.page - 1}" ${ex.page <= 1 ? "disabled" : ""}>${t("sample.prev")}</button>
        <span>${t("sample.tubesInfo", { a: fmtInt(a), b: fmtInt(b), n: fmtInt(x.n_rows) })}</span>
        <button type="button" data-sub="${x.material_id}" data-page="${ex.page + 1}" ${ex.page >= pages ? "disabled" : ""}>${t("sample.next")}</button>
      </div></div>`;
}

async function loadSystemTubes(mid, page) {
  const f = state.sampleFilter || readFilter();
  const ex = state.expanded.get(mid) || { page: 1, rows: [] };
  ex.loading = true;
  state.expanded.set(mid, ex);
  renderSystems();
  try {
    const data = await post("/sample", f, `?limit=${SUB}&offset=${(page - 1) * SUB}&material_id=${mid}`);
    if (!state.expanded.has(mid)) return;
    Object.assign(ex, { page, rows: data.tubes || [], loading: false });
  } catch (e) {
    ex.loading = false;
    toast(e.message, "err");
  }
  renderSystems();
}

function onSampleClick(ev) {
  const tog = ev.target.closest(".cat-toggle");
  if (tog) {
    const mid = +tog.dataset.mid;
    if (state.expanded.has(mid)) { state.expanded.delete(mid); renderSystems(); }
    else loadSystemTubes(mid, 1);
    return;
  }
  const pg = ev.target.closest("button[data-sub]");
  if (pg && !pg.disabled) loadSystemTubes(+pg.dataset.sub, +pg.dataset.page);
}

/* ── Tabela periódica ───────────────────────────────────────────────────── */
const SYMBOLS = ("H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr " +
  "Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In " +
  "Sn Sb Te I Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re " +
  "Os Ir Pt Au Hg Tl Pb Bi Po At Rn Fr Ra Ac Th Pa U Np Pu Am Cm Bk Cf Es Fm " +
  "Md No Lr Rf Db Sg Bh Hs Mt Ds Rg Cn Nh Fl Mc Lv Ts Og").split(" ");
const SYMBOL_BY_LOWER = Object.fromEntries(SYMBOLS.map((s) => [s.toLowerCase(), s]));

function ptPosition(z) {
  if (z === 1) return [1, 1];
  if (z === 2) return [1, 18];
  const short = (i) => (i < 2 ? i + 1 : i + 11);
  if (z <= 10) return [2, short(z - 3)];
  if (z <= 18) return [3, short(z - 11)];
  if (z <= 36) return [4, z - 18];
  if (z <= 54) return [5, z - 36];
  const long = (start, row, frow) => {
    const i = z - start;
    if (i < 2) return [row, i + 1];
    if (i <= 16) return [frow, i + 1];
    return [row, i - 13];
  };
  return z <= 86 ? long(55, 6, 9) : long(87, 7, 10);
}

const split = (s) => s.split(" ");
const CLASSES = [
  // "Todos" exclui tudo o que não está em Contém.  Para ficar só com carbono,
  // "Todos" e depois um clique no C, que sai da exclusão.
  { key: "pt.cAll", els: SYMBOLS },
  { key: "pt.cTM", els: split("Sc Ti V Cr Mn Fe Co Ni Cu Zn Y Zr Nb Mo Tc Ru Rh Pd Ag Cd Hf Ta W Re Os Ir Pt Au Hg") },
  { key: "pt.cLan", els: split("La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu") },
  { key: "pt.cAct", els: split("Ac Th Pa U Np Pu Am Cm Bk Cf Es Fm Md No Lr") },
  { key: "pt.cAlk", els: split("Li Na K Rb Cs Fr") },
  { key: "pt.cAlkE", els: split("Be Mg Ca Sr Ba Ra") },
  { key: "pt.cHal", els: split("F Cl Br I At") },
];
const classPresent = (c) => c.els.filter((e) => state.elementCounts[e]);

/** Por que um elemento não pode ir para o campo `mode`, ou null se pode. */
function blockReason(sym, mode) {
  const total = state.elementCounts[sym];
  if (!total) return "absent";
  if (elState.has(sym)) return null;
  const live = state.liveCounts;
  const n = live ? (live[sym] || 0) : total;
  if (n === 0) return "dim";
  if (mode === "out" && state.liveTotal != null && n >= state.liveTotal) return "all";
  return null;
}

function setElement(sym, mode) {
  if (elState.get(sym) === mode) { elState.delete(sym); return true; }
  const why = blockReason(sym, mode);
  if (why) { warn(why, sym); return false; }
  elState.set(sym, mode);
  return true;
}

let warnTimer = null;
function warn(why, sym) {
  const help = $("pt-help");
  help.textContent = t({ absent: "pt.warnAbsent", dim: "pt.warnDim", all: "pt.warnAll", unknown: "pt.warnUnknown" }[why], { s: sym });
  help.classList.add("is-warn");
  clearTimeout(warnTimer);
  warnTimer = setTimeout(() => { help.textContent = t("pt.help"); help.classList.remove("is-warn"); }, 4000);
}

function setMode(mode) {
  state.mode = mode;
  document.querySelectorAll(".pt-field").forEach((f) => f.classList.toggle("is-active", f.dataset.mode === mode));
  refreshPTable();
}

/** Liga ou desliga a exclusão.  Desligar tira os elementos excluídos. */
function setExclude(on) {
  $("pt-exclude-toggle").checked = on;
  $("pt-field-out").classList.toggle("hidden", !on);
  if (on) {
    setMode("out");
  } else {
    [...elState].forEach(([k, v]) => { if (v === "out") elState.delete(k); });
    setMode("in");
    elementsChanged();
  }
}

function elementsChanged() { renderTokens(); refreshPTable(); scheduleCount(); }

/** Excluir uma classe marca os elementos dela como excluídos (fichas). */
function classOutOn(c) {
  const present = c.els.filter((e) => state.elementCounts[e] && elState.get(e) !== "in");
  return present.length > 0 && present.every((e) => elState.get(e) === "out");
}

function toggleClassIn(c) {
  if (state.classIn.has(c.key)) {
    state.classIn.delete(c.key);
  } else {
    state.classIn.add(c.key);
    classPresent(c).forEach((e) => { if (elState.get(e) === "out") elState.delete(e); });
  }
  elementsChanged();
}

function toggleClassOut(c) {
  const present = classPresent(c);
  if (classOutOn(c)) {
    present.forEach((e) => { if (elState.get(e) === "out") elState.delete(e); });
  } else {
    state.classIn.delete(c.key);
    if (!$("pt-exclude-toggle").checked) setExclude(true);
    present.forEach((e) => { if (!elState.has(e)) elState.set(e, "out"); });
  }
  elementsChanged();
}

function buildPTable() {
  const grid = $("ptable");
  grid.innerHTML = "";
  for (const [row, col, label] of [[6, 3, "57–71"], [7, 3, "89–103"]]) {
    const g = document.createElement("div");
    g.className = "pt-gap";
    g.style.gridRow = row;
    g.style.gridColumn = col;
    g.textContent = label;
    grid.appendChild(g);
  }
  SYMBOLS.forEach((sym, k) => {
    const [row, col] = ptPosition(k + 1);
    const b = document.createElement("button");
    b.type = "button";
    b.className = "pt-el";
    b.dataset.sym = sym;
    b.dataset.z = k + 1;
    b.style.gridRow = row;
    b.style.gridColumn = col;
    b.innerHTML = `<span class="pt-el__z">${k + 1}</span>${sym}`;
    if (!state.elementCounts[sym]) b.disabled = true;
    grid.appendChild(b);
  });
  grid.addEventListener("click", (ev) => {
    const b = ev.target.closest(".pt-el");
    if (!b || b.disabled) return;
    if (setElement(b.dataset.sym, state.mode)) elementsChanged();
  });

  [["pt-classes-out", toggleClassOut]].forEach(([id, fn]) => {
    const box = $(id);
    CLASSES.forEach((c) => {
      const b = document.createElement("button");
      b.type = "button";
      b.className = "pill";
      b.dataset.key = c.key;
      b.addEventListener("click", () => fn(c));
      box.appendChild(b);
    });
  });
}

function refreshPTable() {
  const live = state.liveCounts;
  document.querySelectorAll("#ptable .pt-el").forEach((b) => {
    const sym = b.dataset.sym;
    const st = elState.get(sym);
    b.classList.toggle("is-in", st === "in");
    b.classList.toggle("is-out", st === "out");
    b.setAttribute("aria-pressed", st ? "true" : "false");
    const total = state.elementCounts[sym];
    const why = total ? blockReason(sym, state.mode) : "absent";
    b.classList.toggle("is-dim", Boolean(total) && !st && Boolean(why));
    const head = `${sym} (${b.dataset.z})`;
    if (!total) b.title = `${head} · ${t("pt.absent")}`;
    else if (why === "dim") b.title = `${head} · ${t("pt.dim")}`;
    else if (why === "all") b.title = `${head} · ${t("pt.warnAll", { s: sym })}`;
    else if (live) b.title = `${head} · ${fmtInt(live[sym] || 0)} ${t("pt.withSel")}`;
    else b.title = `${head} · ${fmtInt(total)} ${t("pt.materials")}`;
  });
  CLASSES.forEach((c) => {
    const bin = document.querySelector(`#pt-classes-in [data-key="${c.key}"]`);
    const bout = document.querySelector(`#pt-classes-out [data-key="${c.key}"]`);
    if (bin) { bin.textContent = t(c.key); bin.classList.toggle("is-on", state.classIn.has(c.key)); bin.setAttribute("aria-pressed", String(state.classIn.has(c.key))); }
    if (bout) { bout.textContent = t(c.key); bout.classList.toggle("is-on", classOutOn(c)); bout.setAttribute("aria-pressed", String(classOutOn(c))); }
  });
}

/* Campos digitáveis: fichas para o que foi aceito e uma caixa de texto no fim.
 * Espaço, vírgula, Enter ou sair do campo aceitam; Backspace vazio tira a
 * última ficha. */
function renderTokens() {
  ["in", "out"].forEach((mode) => {
    const box = $(`pt-box-${mode}`);
    const input = $(`pt-input-${mode}`);
    box.querySelectorAll(".tok").forEach((x) => x.remove());
    const syms = [...elState].filter(([, k]) => k === mode).map(([sym]) => sym);
    // Muitos excluídos viram uma ficha só ("Todos" deixaria 85 fichas no campo).
    if (mode === "out" && syms.length > 10) {
      const tok = document.createElement("span");
      tok.className = "tok tok--out";
      tok.title = syms.join(" ");
      tok.innerHTML = `${t("pt.outMany", { n: syms.length })}<button type="button" aria-label="${t("pt.clear")}" data-clear-out="1">×</button>`;
      box.insertBefore(tok, input);
      return;
    }
    syms.map((sym) => [sym]).forEach(([sym]) => {
      const tok = document.createElement("span");
      tok.className = `tok tok--${mode}`;
      tok.innerHTML = `${sym}<button type="button" aria-label="${sym}" data-sym="${sym}">×</button>`;
      box.insertBefore(tok, input);
    });
  });
}

function commitTyped(mode) {
  const input = $(`pt-input-${mode}`);
  const parts = input.value.split(/[\s,;]+/).filter(Boolean);
  input.value = "";
  let changed = false;
  for (const raw of parts) {
    const sym = SYMBOL_BY_LOWER[raw.toLowerCase()];
    if (!sym) { warn("unknown", raw); continue; }
    if (elState.get(sym) === mode) continue;
    if (setElement(sym, mode)) changed = true;
  }
  if (changed) elementsChanged();
}

function wireTokenFields() {
  ["in", "out"].forEach((mode) => {
    const field = $(`pt-field-${mode}`), box = $(`pt-box-${mode}`), input = $(`pt-input-${mode}`);
    field.addEventListener("mousedown", () => setMode(mode));
    box.addEventListener("click", (ev) => {
      if (ev.target.closest("button[data-clear-out]")) {
        [...elState].forEach(([k, v]) => { if (v === "out") elState.delete(k); });
        elementsChanged();
        return;
      }
      const x = ev.target.closest("button[data-sym]");
      if (x) { elState.delete(x.dataset.sym); elementsChanged(); return; }
      input.focus();
    });
    input.addEventListener("focus", () => setMode(mode));
    input.addEventListener("keydown", (ev) => {
      if (ev.key === "Enter" || ev.key === "," || ev.key === " ") { ev.preventDefault(); commitTyped(mode); }
      else if (ev.key === "Backspace" && !input.value) {
        const last = [...elState].filter(([, k]) => k === mode).pop();
        if (last) { elState.delete(last[0]); elementsChanged(); }
      }
    });
    input.addEventListener("input", () => { if (/[\s,;]/.test(input.value)) commitTyped(mode); });
    input.addEventListener("blur", () => commitTyped(mode));
  });
  $("pt-exclude-toggle").addEventListener("change", (ev) => setExclude(ev.target.checked));
}

let elSeq = 0;
async function refreshElementCounts() {
  // Sempre pergunta: os filtros de tubo (ligações espúrias, diâmetro, átomos,
  // cavidade, quiralidade) também decidem quais elementos ainda dão tubos.
  const f = readFilter();
  const seq = (elSeq += 1);
  if (filterIsEmpty()) {
    state.liveCounts = null;
    state.liveTotal = null;
    refreshPTable();
    return;
  }
  try {
    const r = await post("/elements", f);
    if (seq !== elSeq) return;
    state.liveCounts = r.counts;
    state.liveTotal = r.n_materials;
  } catch (_) {
    if (seq !== elSeq) return;
    state.liveCounts = null;
    state.liveTotal = null;
  }
  refreshPTable();
  renderPager();        // na amostra agrupada, o total de páginas vem daqui
}

/* ── Figura das classes de rede ─────────────────────────────────────────── */
function openFigure() {
  const img = $("fig-img");
  if (!img.getAttribute("src")) img.setAttribute("src", "static/img/fig_bravais_to_tube.png");
  img.alt = t("fig.title");
  $("fig-overlay").classList.add("open");
}
function closeFigure() { $("fig-overlay").classList.remove("open"); }

/* ── Pedido, montagem na hora e acompanhamento ──────────────────────────── */
async function sendRequest() {
  const c = state.count;
  if (!c || !c.n_tubes) { toast(t("msg.empty")); return; }
  const email = ($("req-email").value || "").trim();
  if (!/^[^@\s]+@[^@\s]+\.[^@\s]+$/.test(email)) { toast(t("msg.needEmail"), "err"); return; }
  const btn = $("req-send");
  btn.disabled = true;
  try {
    const data = await post("/request", { filter: readFilter(), email, name: ($("req-name").value || "").trim() });
    if (data.mode === "instant") {
      toast(t("msg.instant"), "ok");
      openBuild(data.protocol);
      return;
    }
    if (data.cached) toast(t("msg.cached"), "ok");
    else if (data.mail_sent === false) toast(t("msg.mailFail"), "info");
    else toast(t("msg.sent"), "ok");
    startTracking(data.protocol);
  } catch (e) {
    toast(e.message, "err");
  } finally {
    btn.disabled = false;
  }
}

function openBuild(proto) {
  state.build = { proto, downloaded: false };
  $("build-proto").textContent = proto;
  $("build-fill").style.width = "0%";
  $("build-state").textContent = t("build.wait");
  $("build-dl").classList.add("hidden");
  $("build-overlay").classList.add("open");
  $("track-proto").value = proto;
  pollBuild();
}

function closeBuild() {
  clearTimeout(state.buildTimer);
  state.build = null;
  $("build-overlay").classList.remove("open");
  if ($("track-proto").value) pollTrack();
}

async function pollBuild() {
  clearTimeout(state.buildTimer);
  const b = state.build;
  if (!b) return;
  try {
    const r = await fetch(`${API}/status/${encodeURIComponent(b.proto)}`);
    const d = await r.json();
    if (state.build !== b) return;
    if (d.state === "done") {
      $("build-fill").style.width = "100%";
      $("build-state").textContent = t("build.done");
      const a = $("build-dl");
      a.href = `${API}/download/${encodeURIComponent(b.proto)}`;
      a.classList.remove("hidden");
      if (!b.downloaded) { b.downloaded = true; window.location.href = a.href; }
      return;
    }
    if (d.state === "failed" || d.state === "expired") {
      $("build-state").textContent = t("build.failed", { m: d.message || "" });
      return;
    }
    $("build-fill").style.width = `${d.percent || 0}%`;
    $("build-state").textContent = d.state === "running" && d.progress
      ? t("build.running", { d: fmtInt(d.progress), n: fmtInt(d.n_tubes) })
      : t("build.wait");
  } catch (e) {
    $("build-state").textContent = e.message;
  }
  state.buildTimer = setTimeout(pollBuild, 1000);
}

function startTracking(protocol) {
  $("track-proto").value = protocol;
  $("cat-track").scrollIntoView({ behavior: "smooth", block: "center" });
  pollTrack();
}

async function pollTrack() {
  clearTimeout(state.trackTimer);
  const proto = ($("track-proto").value || "").trim();
  if (!proto) return;
  const out = $("track-out");
  try {
    const r = await fetch(`${API}/status/${encodeURIComponent(proto)}`);
    if (r.status === 404) { out.innerHTML = `<p class="cat-state">${t("msg.noProto")}</p>`; return; }
    const d = await r.json();
    const pct = d.percent || 0;
    const label = d.state === "queued" && d.mode === "instant" ? t("state.localQueued") : t(`state.${d.state}`);
    out.innerHTML = `
      <p><span class="cat-proto">${d.protocol}</span> <span class="cat-state">· ${label}</span></p>
      <div class="cat-bar"><div class="cat-bar__fill" style="width:${pct}%"></div></div>
      <p class="cat-state">${fmtInt(d.progress || 0)} / ${fmtInt(d.n_tubes)} · ${fmtNum(pct, 1)} %</p>
      ${d.state === "done"
        ? `<p><a class="btn-tonal" href="${API}/download/${d.protocol}">${t("dl.link")} (${fmtBytes(d.bytes)})</a>
           <span class="cat-state">${t("dl.expires")} ${(d.expires || "").slice(0, 10)}</span></p>`
        : ""}`;
    if (d.state === "queued" || d.state === "running") state.trackTimer = setTimeout(pollTrack, 5000);
  } catch (e) {
    out.innerHTML = `<p class="cat-state">${e.message}</p>`;
  }
}

/* ── Como funciona: balões empacotados ──────────────────────────────────────
 * Objetivo: o menor retângulo que contém todos os balões, sem sobra no fim das
 * colunas nem espaço vazio dentro dos balões.
 *
 * A altura de um balão numa largura w segue h(w) ≈ a/w + b (as linhas de texto
 * caem com 1/w); a e b saem de duas medidas.  Para cada número de colunas k que
 * cabe na largura, e para cada divisão dos balões em k colunas mantendo a ordem
 * de leitura, as larguras que igualam as alturas saem de
 * Σ_c A_c/(H − B_c) = largura disponível, resolvido por bisseção em H.  Fica a
 * divisão de menor H com todas as colunas acima da largura mínima.  O layout só
 * é refeito quando a largura muda ou o idioma troca. */
const EXPLAIN = { width: 0, articles: null };

function* cutsOf(n, k, start = 1) {
  // posições de corte crescentes em 1..n-1, k-1 delas
  if (k === 1) { yield []; return; }
  for (let i = start; i <= n - (k - 1); i += 1) {
    for (const rest of cutsOf(n, k - 1, i + 1)) yield [i, ...rest];
  }
}

/** Todas as divisões dos balões em colunas, em ordem de leitura. */
function allPartitions(n, maxK) {
  const out = [];
  for (let k = 1; k <= maxK; k += 1) for (const cuts of cutsOf(n, k)) out.push([0, ...cuts, n]);
  return out;
}

/** Altura real de cada balão em cada largura, de STEP em STEP px, numa passada
 *  só: todos os clones entram fora da tela de uma vez e o navegador calcula o
 *  layout uma vez (ler a altura depois de cada mudança forçava milhares de
 *  recálculos, 3 s na janela de 1280 px).  A altura de uma coluna é a soma dos
 *  seus balões, e a tabela usa a largura imediatamente abaixo, então nunca
 *  promete menos altura do que a real. */
function measureTable(box, arts, minW, maxW, step) {
  const widths = [];
  for (let w = minW; w < maxW; w += step) widths.push(w);
  widths.push(maxW);
  const probe = document.createElement("div");
  probe.className = "cat-explain cat-explain--probe";
  const cells = [];
  const frag = document.createDocumentFragment();
  arts.forEach((a, i) => widths.forEach((w, j) => {
    const c = a.cloneNode(true);
    c.style.width = `${w}px`;
    c.style.margin = "0";
    frag.appendChild(c);
    cells.push([i, j, c]);
  }));
  probe.appendChild(frag);
  box.parentElement.appendChild(probe);
  const table = arts.map(() => new Float64Array(widths.length));
  cells.forEach(([i, j, c]) => { table[i][j] = c.getBoundingClientRect().height; });
  probe.remove();
  return { widths, table, step, minW };
}

/** Menor altura da caixa em que cada coluna da divisão cabe na largura dada. */
function exactPlan(T, bounds, W, GAP, MIN_W) {
  const k = bounds.length - 1;
  const avail = W - GAP * (k - 1);
  const maxW = avail - MIN_W * (k - 1);
  const idxOf = (w) => Math.max(0, Math.min(T.widths.length - 1, Math.floor((w - T.minW) / T.step)));
  const maxIdx = idxOf(maxW);
  const colH = (c, j) => {
    let sum = GAP * (bounds[c + 1] - bounds[c] - 1);
    for (let i = bounds[c]; i < bounds[c + 1]; i += 1) sum += T.table[i][j];
    return sum;
  };
  const minIdxFor = (c, H) => {
    if (colH(c, maxIdx) > H) return -1;
    let lo = -1, hi = maxIdx;           // colH(hi) <= H
    while (hi - lo > 1) {
      const mid = (lo + hi) >> 1;
      if (colH(c, mid) <= H) hi = mid; else lo = mid;
    }
    return hi;
  };
  const cols = [...Array(k).keys()];
  const sumFor = (H) => {
    let total = 0;
    for (const c of cols) {
      const j = minIdxFor(c, H);
      if (j < 0) return Infinity;
      total += T.widths[j];
    }
    return total;
  };
  let Hlo = Math.max(...cols.map((c) => colH(c, maxIdx))) - 0.5;
  let Hhi = Math.max(...cols.map((c) => colH(c, 0)));
  while (Hhi - Hlo > 0.5) {
    const mid = (Hlo + Hhi) / 2;
    if (sumFor(mid) <= avail) Hhi = mid; else Hlo = mid;
  }
  const idx = cols.map((c) => minIdxFor(c, Hhi));
  const widths = idx.map((j) => T.widths[j]);
  // A sobra de largura vai, passo a passo, para a coluna mais alta.
  let slack = avail - widths.reduce((x, y) => x + y, 0);
  while (slack >= T.step) {
    const tallest = cols.reduce((m, c) => (colH(c, idx[c]) > colH(m, idx[m]) ? c : m), 0);
    if (idx[tallest] >= T.widths.length - 1) break;
    idx[tallest] += 1;
    widths[tallest] += T.step;
    slack -= T.step;
  }
  if (slack > 0) widths[0] += slack;
  const heights = cols.map((c) => colH(c, idx[c]));
  const H = Math.max(...heights);
  return { bounds, widths, H, gaps: heights.reduce((acc, x) => acc + (H - x), 0) };
}

function applyExplain(box, arts, plan) {
  box.innerHTML = "";
  for (let c = 0; c < plan.widths.length; c += 1) {
    const col = document.createElement("div");
    col.className = "cat-explain__col";
    col.style.flex = `${plan.widths[c].toFixed(2)} 1 0`;
    arts.slice(plan.bounds[c], plan.bounds[c + 1]).forEach((a) => col.appendChild(a));
    box.appendChild(col);
  }
}

function layoutExplain(force) {
  const box = $("cat-explain-grid");
  if (!box || !box.clientWidth) return;
  const W = box.clientWidth;
  if (!force && W === EXPLAIN.width) return;
  EXPLAIN.width = W;
  if (!EXPLAIN.articles) EXPLAIN.articles = [...box.querySelectorAll("article")];
  const arts = EXPLAIN.articles;
  const GAP = 12, MIN_W = 230;
  const t0 = performance.now();
  const maxK = Math.max(1, Math.min(arts.length, Math.floor((W + GAP) / (MIN_W + GAP))));
  const T = measureTable(box, arts, MIN_W, W, 4);
  let best = null;
  for (const bounds of allPartitions(arts.length, maxK)) {
    const plan = exactPlan(T, bounds, W, GAP, MIN_W);
    // Retângulo e espaço vazio juntos: a altura da caixa mais a sobra média por
    // coluna, que é o vazio que cada balão da coluna vai carregar.  Só a altura
    // trocava 6 px de caixa por 56 px de vazio (janela de 1100 px).
    plan.score = plan.H + plan.gaps / (bounds.length - 1);
    if (!best || plan.score < best.score) best = plan;
  }
  applyExplain(box, arts, best);
  EXPLAIN.best = best;
  EXPLAIN.ms = Math.round(performance.now() - t0);
}
let explainTimer = null;
const scheduleExplain = (force) => {
  clearTimeout(explainTimer);
  explainTimer = setTimeout(() => layoutExplain(force === true), 60);
};

/* ── Idioma ─────────────────────────────────────────────────────────────── */
function applyLang(lang) {
  LANG = I18N[lang] ? lang : "pt";
  document.documentElement.lang = LANG;
  document.querySelectorAll("[data-i18n]").forEach((el) => {
    const k = el.getAttribute("data-i18n");
    if (I18N[LANG][k]) el.textContent = I18N[LANG][k];
  });
  $("track-proto").placeholder = t("track.placeholder");
  // O style.css do construtor pinta o idioma ativo pela classe .active.
  document.querySelectorAll(".lang-toggle button").forEach((b) => {
    const on = b.dataset.lang === LANG;
    b.classList.toggle("active", on);
    b.setAttribute("aria-pressed", String(on));
  });
  try { localStorage.setItem("ntbuilder-lang", LANG); } catch (_) {}
  scheduleExplain(true);
  if (!state.meta) return;
  document.querySelectorAll("#f-lattices [data-v]").forEach((b) => { b.textContent = t(`lat.${b.dataset.v}`); });
  renderStats();
  sliderReadouts();
  refreshPTable();
  renderCount();
  renderCurrent();
}

/* ── Montagem ───────────────────────────────────────────────────────────── */
async function boot() {
  try { document.documentElement.setAttribute("data-theme", localStorage.getItem("ntbuilder-theme") || "light"); } catch (_) {}
  $("theme-toggle-btn").addEventListener("click", () => {
    const cur = document.documentElement.getAttribute("data-theme") === "dark" ? "light" : "dark";
    document.documentElement.setAttribute("data-theme", cur);
    try { localStorage.setItem("ntbuilder-theme", cur); } catch (_) {}
  });
  let lang = "pt";
  try { lang = localStorage.getItem("ntbuilder-lang") || "pt"; } catch (_) {}
  applyLang(lang);
  document.querySelectorAll(".lang-toggle button").forEach((b) => b.addEventListener("click", () => applyLang(b.dataset.lang)));
  // Só a largura importa: o próprio layout muda a altura da caixa, e observar
  // qualquer mudança de tamanho faria o layout se refazer sem parar.
  // A fonte da página chega depois: o texto muda de altura sem a largura mudar.
  if (document.fonts && document.fonts.ready) document.fonts.ready.then(() => scheduleExplain(true));
  if (window.ResizeObserver) new ResizeObserver(() => scheduleExplain()).observe($("cat-explain-grid"));
  else window.addEventListener("resize", () => scheduleExplain());

  let meta;
  try { meta = await (await fetch(`${API}/meta`)).json(); } catch (e) { toast(e.message, "err"); return; }
  if (meta.detail) { toast(meta.detail, "err"); return; }
  state.meta = meta;
  state.elementCounts = meta.element_counts || {};
  // Recursos que a API anuncia.  Uma API anterior ignora classes e páginas, e
  // aí os controles delas somem em vez de mostrar números errados.
  const feats = new Set(meta.features || []);
  state.grouped = feats.has("systems");
  if (!feats.has("pages")) document.querySelector(".cat-pager").classList.add("hidden");
  state.canOpen = feats.has("layer");
  if (meta.instant_max_seconds) state.instantMax = meta.instant_max_seconds;

  // O diâmetro vai até o maior tubo que existe no banco.
  const lim = meta.limits || {};
  state.dMax = Math.max(30, Math.ceil(lim.d_max || 80));
  ["f-d-min", "f-d-max"].forEach((id) => { $(id).max = state.dMax; });
  $("f-d-max").value = state.dMax;

  const dbs = $("f-sources");
  const order = (meta.by_source_stats || []).slice().sort((a, b) => b.n_materials - a.n_materials).map((s) => s.source);
  (order.length ? order : meta.by_source.map((s) => s.source)).forEach((source) => {
    const b = document.createElement("button");
    b.type = "button";
    b.className = "db-toggle is-on";
    b.dataset.v = source;
    b.setAttribute("aria-pressed", "true");
    b.innerHTML = `<span class="db-toggle__check">✓</span><span class="db-toggle__body"><span class="db-toggle__name"></span><span class="db-toggle__n"></span></span>`;
    dbs.appendChild(b);
  });
  const lats = $("f-lattices");
  meta.by_lattice.forEach((s) => {
    const b = document.createElement("button");
    b.type = "button";
    b.className = "pill is-on";
    b.dataset.v = s.lattice_type;
    lats.appendChild(b);
  });

  buildPTable();
  wireTokenFields();
  wirePager();
  $("sample-table").addEventListener("click", onSampleClick);
  wirePills("f-sources", () => { renderStats(); scheduleCount(); });
  wirePills("f-lattices", scheduleCount);
  wirePills("f-chirality", scheduleCount);
  wirePills("f-formats", () => { renderSummary(); renderCurrent(); scheduleCount(); });
  wireSeg("f-clean", scheduleCount);
  wireDual("f-d-min", "f-d-max", scheduleCount);
  wireDual("f-t-min", "f-t-max", scheduleCount);
  ["f-din", "f-atoms", "f-eps"].forEach((id) => $(id).addEventListener("input", () => { sliderReadouts(); scheduleCount(); }));
  ["f-vacuum", "f-per-mat"].forEach((id) => $(id).addEventListener("input", renderSummary));
  $("pt-clear").addEventListener("click", () => { elState.clear(); state.classIn.clear(); elementsChanged(); });
  $("btn-reset").addEventListener("click", () => location.reload());
  $("req-send").addEventListener("click", sendRequest);
  $("btn-track").addEventListener("click", pollTrack);
  $("lat-help").addEventListener("click", openFigure);
  $("fig-close").addEventListener("click", closeFigure);
  $("fig-overlay").addEventListener("click", (ev) => { if (ev.target.id === "fig-overlay") closeFigure(); });
  $("build-close").addEventListener("click", closeBuild);
  document.addEventListener("keydown", (ev) => { if (ev.key === "Escape") closeFigure(); });

  applyLang(LANG);
  renderTokens();
  refreshCount();

  // O link do e-mail chega como /catalogo?protocolo=NTB-...
  const q = new URLSearchParams(location.search);
  const proto = q.get("protocolo") || q.get("protocol");
  if (proto) startTracking(proto);
}

boot();
